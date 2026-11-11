//! One-pass slicing of a burn into a regular tile grid.
//!
//! [`slice_tiles`] routes every record of a whole-grid burn to the tile
//! that holds it, in a single pass: a record in cell `(row, col)` goes to
//! tile `(row / tile_nrow, col / tile_ncol)` with local indices
//! `(row % tile_nrow, col % tile_ncol)`, and a run is split wherever it
//! crosses a tile column seam. Nothing is recomputed, so each tile equals
//! [`crate::crop`] of the full result to that tile's window, record for
//! record, and the materialized tiles put back together equal the
//! materialized whole grid bit for bit.
//!
//! Cost is O(records + tiles). Calling [`crate::crop`] once per tile gives
//! the same tables but scans every record for every tile, which is
//! O(records × tiles); use `crop` for one window and `slice_tiles` for a
//! tile grid.
//!
//! [`slice_tiles_par`] gives the same tiles using several threads (std
//! only): records are partitioned by band of tile rows, keeping their
//! order, and the bands are sliced in parallel.

use std::ops::Range;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::grid::GridSpec;
use crate::output::{BurnResult, GridEdge, GridLine, GridPoint, GridRun};

/// One non-empty tile of a [`slice_tiles`] call.
#[derive(Clone, Debug, PartialEq)]
pub struct Tile {
    /// Tile row in the tile grid (0 at the top).
    pub tile_row: u32,
    /// Tile column in the tile grid.
    pub tile_col: u32,
    /// Row of the tile's first cell in the whole grid.
    pub row_off: u32,
    /// Column of the tile's first cell in the whole grid.
    pub col_off: u32,
    /// The tile's own grid. Tiles on the right and bottom edges are
    /// smaller when the tile shape does not divide the grid.
    pub grid: GridSpec,
    /// The tile's records, indexed on `grid`, in burn order. `notes` is
    /// left empty: notes are per geometry, so read them from the whole
    /// result.
    pub result: BurnResult,
}

/// Slice `result` (burned into `grid`) into tiles of `tile_ncol` ×
/// `tile_nrow` cells, anchored at the grid's top-left cell.
///
/// Returns only the tiles that hold at least one record, in row-major
/// tile order. The tile grid has `ceil(ncol / tile_ncol)` columns and
/// `ceil(nrow / tile_nrow)` rows.
///
/// # Panics
///
/// If `tile_ncol` or `tile_nrow` is zero, or if a record lies outside
/// `grid` (which `burn` never produces for the same `grid`).
///
/// # Examples
///
/// ```
/// use controlledburn::{burn, slice_tiles, BurnOptions, Coord, Geometry,
///     GridSpec, Polygon};
///
/// let square = Geometry::Polygon(Polygon::new(vec![vec![
///     Coord::new(2.5, 4.5), Coord::new(6.5, 4.5), Coord::new(6.5, 8.5),
///     Coord::new(2.5, 8.5), Coord::new(2.5, 4.5),
/// ]]));
/// let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
/// let r = burn(&[square], &grid, BurnOptions::default()).unwrap();
///
/// let tiles = slice_tiles(&r, &grid, 4, 4);
/// // the square touches cell rows 1..=5 and columns 2..=6, so tile rows
/// // 0..=1 and tile columns 0..=1
/// assert_eq!(tiles.len(), 4);
/// assert_eq!((tiles[0].grid.ncol, tiles[0].grid.nrow), (4, 4));
/// ```
pub fn slice_tiles(result: &BurnResult, grid: &GridSpec, tile_ncol: u32, tile_nrow: u32) -> Vec<Tile> {
    let s = Shape::new(grid, tile_ncol, tile_nrow);
    s.slice(result, 0..s.ntr, None)
}

/// Below this many records `slice_tiles_par` slices on the calling thread:
/// starting threads costs more than the slice.
const PAR_MIN_RECORDS: usize = 1 << 15;

/// Bands of tile rows per thread: enough that one dense band does not
/// leave the other threads idle.
const BANDS_PER_THREAD: usize = 8;

/// [`slice_tiles`] on `threads` threads (0: all available). The tiles are
/// the same, in the same order, with the same records in the same order.
///
/// The records are split into one contiguous chunk per thread; each chunk
/// counts its records per band of tile rows, then writes their indices
/// into its place in each band's list, so a band's list is in burn order;
/// then the bands are sliced in parallel, each into tile rows no other
/// band touches. The extra memory is one `u32` per record.
///
/// With one thread, a single tile row, or fewer than about 32,000 records
/// it runs [`slice_tiles`] on the calling thread. Threads come from
/// `std::thread::scope`; on a target without threads (such as
/// `wasm32-unknown-unknown`) pass `threads = 1`, or use [`slice_tiles`].
///
/// # Panics
///
/// As [`slice_tiles`].
pub fn slice_tiles_par(
    result: &BurnResult,
    grid: &GridSpec,
    tile_ncol: u32,
    tile_nrow: u32,
    threads: usize,
) -> Vec<Tile> {
    slice_tiles_par_from(result, grid, tile_ncol, tile_nrow, threads, PAR_MIN_RECORDS)
}

/// `slice_tiles_par`, going parallel from `min_records` records (the
/// tests force the parallel path on small inputs).
fn slice_tiles_par_from(
    result: &BurnResult,
    grid: &GridSpec,
    tile_ncol: u32,
    tile_nrow: u32,
    threads: usize,
    min_records: usize,
) -> Vec<Tile> {
    let s = Shape::new(grid, tile_ncol, tile_nrow);
    let threads = if threads == 0 { std::thread::available_parallelism().map_or(1, |n| n.get()) } else { threads };
    let lens = [result.runs.len(), result.edges.len(), result.lines.len(), result.points.len()];
    let total: usize = lens.iter().sum();
    if threads < 2 || s.ntr < 2 || total < min_records || lens.iter().any(|&n| n > u32::MAX as usize) {
        return s.slice(result, 0..s.ntr, None);
    }

    // Band b holds tile rows band_start(b)..band_start(b + 1); a tile row's
    // band comes from a table, so a record costs one division.
    let ntr = s.ntr;
    let nbands = (threads * BANDS_PER_THREAD).min(ntr).min(u16::MAX as usize);
    let band_start = |b: usize| (b * ntr + nbands - 1) / nbands;
    let band_of_tile_row: Vec<u16> = (0..ntr).map(|t| (t * nbands / ntr) as u16).collect();
    let nchunks = threads;
    let chunk = |len: usize, c: usize| (c * len / nchunks)..((c + 1) * len / nchunks);

    // 1. Each record's band (checking the record is in the grid), and per
    //    chunk the records per band.
    let mut bands: [Vec<u16>; 4] = std::array::from_fn(|t| vec![0u16; lens[t]]);
    let counts: Vec<[Vec<usize>; 4]> = {
        let (s, band_of_tile_row) = (&s, &band_of_tile_row);
        let mut pieces: Vec<[&mut [u16]; 4]> = (0..nchunks).map(|_| Default::default()).collect();
        for (t, v) in bands.iter_mut().enumerate() {
            let mut rest: &mut [u16] = v;
            for (c, p) in pieces.iter_mut().enumerate() {
                let (head, tail) = std::mem::take(&mut rest).split_at_mut(chunk(lens[t], c).len());
                p[t] = head;
                rest = tail;
            }
        }
        std::thread::scope(|sc| {
            let hs: Vec<_> = pieces
                .into_iter()
                .enumerate()
                .map(|(c, [b0, b1, b2, b3])| {
                    sc.spawn(move || {
                        let place = |row: i32, col: i32| {
                            s.check(row, col);
                            band_of_tile_row[row as usize / tile_nrow as usize]
                        };
                        [
                            band_chunk(&result.runs[chunk(lens[0], c)], b0, nbands, |r| place(r.row, r.col_start)),
                            band_chunk(&result.edges[chunk(lens[1], c)], b1, nbands, |e| place(e.row, e.col)),
                            band_chunk(&result.lines[chunk(lens[2], c)], b2, nbands, |l| place(l.row, l.col)),
                            band_chunk(&result.points[chunk(lens[3], c)], b3, nbands, |p| place(p.row, p.col)),
                        ]
                    })
                })
                .collect();
            hs.into_iter().map(join).collect()
        })
    };

    // 2. Each table's record indices, band by band and within a band chunk
    //    by chunk, so in burn order. Each chunk writes its own disjoint
    //    pieces.
    let mut idx: [Vec<u32>; 4] = std::array::from_fn(|t| vec![0u32; lens[t]]);
    let mut offsets: [Vec<usize>; 4] = std::array::from_fn(|_| Vec::with_capacity(nbands + 1));
    {
        let mut pieces: Vec<[Vec<&mut [u32]>; 4]> =
            (0..nchunks).map(|_| std::array::from_fn(|_| Vec::with_capacity(nbands))).collect();
        for (t, v) in idx.iter_mut().enumerate() {
            let mut rest: &mut [u32] = v;
            let mut at = 0;
            for b in 0..nbands {
                offsets[t].push(at);
                for (c, cn) in counts.iter().enumerate() {
                    let (head, tail) = std::mem::take(&mut rest).split_at_mut(cn[t][b]);
                    pieces[c][t].push(head);
                    rest = tail;
                    at += cn[t][b];
                }
            }
            offsets[t].push(at);
        }
        let bands = &bands;
        std::thread::scope(|sc| {
            let hs: Vec<_> = pieces
                .into_iter()
                .enumerate()
                .map(|(c, mut mine)| {
                    sc.spawn(move || {
                        for (t, out) in mine.iter_mut().enumerate() {
                            let r = chunk(lens[t], c);
                            let mut pos = vec![0usize; nbands];
                            for (i, &b) in r.clone().zip(&bands[t][r]) {
                                let b = b as usize;
                                out[b][pos[b]] = i as u32;
                                pos[b] += 1;
                            }
                        }
                    })
                })
                .collect();
            hs.into_iter().for_each(join);
        });
    }
    let s = &s;

    // 3. Slice the bands in parallel, then put them back in band order.
    let next = AtomicUsize::new(0);
    let (idx, offsets, next) = (&idx, &offsets, &next);
    let mut parts: Vec<(usize, Vec<Tile>)> = std::thread::scope(|sc| {
        let hs: Vec<_> = (0..threads.min(nbands))
            .map(|_| {
                sc.spawn(move || {
                    let mut mine = Vec::new();
                    loop {
                        let b = next.fetch_add(1, Ordering::Relaxed);
                        if b >= nbands {
                            return mine;
                        }
                        let band = |t: usize| &idx[t][offsets[t][b]..offsets[t][b + 1]];
                        let sel = Sel { runs: band(0), edges: band(1), lines: band(2), points: band(3) };
                        mine.push((b, s.slice(result, band_start(b)..band_start(b + 1), Some(sel))));
                    }
                })
            })
            .collect();
        hs.into_iter().flat_map(join).collect()
    });
    parts.sort_unstable_by_key(|p| p.0);
    parts.into_iter().flat_map(|p| p.1).collect()
}

/// Join a scoped thread, re-raising its panic with its own payload (so
/// the message is the one `slice_tiles` would give).
fn join<T>(h: std::thread::ScopedJoinHandle<'_, T>) -> T {
    h.join().unwrap_or_else(|e| std::panic::resume_unwind(e))
}

/// The band of each record of a chunk (`place`: its band, after checking
/// it lies in the grid), written to `out`; returns the records per band.
fn band_chunk<T>(v: &[T], out: &mut [u16], nbands: usize, place: impl Fn(&T) -> u16) -> Vec<usize> {
    let mut n = vec![0usize; nbands];
    for (x, b) in v.iter().zip(out.iter_mut()) {
        *b = place(x);
        n[*b as usize] += 1;
    }
    n
}

/// The records of one band, as indices into each table, in burn order.
#[derive(Clone, Copy)]
struct Sel<'a> {
    runs: &'a [u32],
    edges: &'a [u32],
    lines: &'a [u32],
    points: &'a [u32],
}

/// Every record of `v`, or the selected ones, in order.
fn each<T>(v: &[T], sel: Option<&[u32]>, mut f: impl FnMut(&T)) {
    match sel {
        None => v.iter().for_each(f),
        Some(ix) => ix.iter().for_each(|&i| f(&v[i as usize])),
    }
}

/// The tile grid over `grid`.
struct Shape<'a> {
    grid: &'a GridSpec,
    tile_ncol: u32,
    tile_nrow: u32,
    /// tile columns and rows
    ntc: usize,
    ntr: usize,
}

impl<'a> Shape<'a> {
    fn new(grid: &'a GridSpec, tile_ncol: u32, tile_nrow: u32) -> Self {
        assert!(tile_ncol > 0 && tile_nrow > 0, "tile shape must be at least 1 x 1");
        let ntc = (grid.ncol as usize + tile_ncol as usize - 1) / tile_ncol as usize;
        let ntr = (grid.nrow as usize + tile_nrow as usize - 1) / tile_nrow as usize;
        Shape { grid, tile_ncol, tile_nrow, ntc, ntr }
    }

    fn check(&self, row: i32, col: i32) {
        let g = self.grid;
        assert!(
            row >= 0 && col >= 0 && (row as u32) < g.nrow && (col as u32) < g.ncol,
            "record ({row}, {col}) is outside the {} x {} grid",
            g.ncol,
            g.nrow
        );
    }

    /// The tiles in tile rows `rows`, from every record (`sel` None) or
    /// the selected ones, all of which lie in those tile rows.
    fn slice(&self, result: &BurnResult, rows: Range<usize>, sel: Option<Sel>) -> Vec<Tile> {
        let grid = self.grid;
        let (tile_ncol, tile_nrow) = (self.tile_ncol, self.tile_nrow);
        let (tc, tr) = (tile_ncol as i32, tile_nrow as i32);
        let ntc = self.ntc;
        let r0 = rows.start;
        let at = |row: i32, col: i32| -> usize {
            self.check(row, col);
            ((row / tr) as usize - r0) * ntc + (col / tc) as usize
        };

        let mut out: Vec<BurnResult> = vec![BurnResult::default(); ntc * rows.len()];
        each(&result.runs, sel.map(|s| s.runs), |run| {
            let lr = run.row % tr;
            let mut c = run.col_start;
            while c < run.col_end {
                let c0 = (c / tc) * tc;
                let ce = (c0 + tc).min(run.col_end);
                out[at(run.row, c)].runs.push(GridRun { row: lr, col_start: c - c0, col_end: ce - c0, id: run.id });
                c = ce;
            }
        });
        each(&result.edges, sel.map(|s| s.edges), |e| {
            out[at(e.row, e.col)].edges.push(GridEdge { row: e.row % tr, col: e.col % tc, ..*e });
        });
        each(&result.lines, sel.map(|s| s.lines), |l| {
            out[at(l.row, l.col)].lines.push(GridLine { row: l.row % tr, col: l.col % tc, ..*l });
        });
        each(&result.points, sel.map(|s| s.points), |p| {
            out[at(p.row, p.col)].points.push(GridPoint { row: p.row % tr, col: p.col % tc, ..*p });
        });

        let (dx, dy) = (grid.dx(), grid.dy());
        out.into_iter()
            .enumerate()
            .filter(|(_, r)| !(r.runs.is_empty() && r.edges.is_empty() && r.lines.is_empty() && r.points.is_empty()))
            .map(|(k, result)| {
                let (ti, tj) = ((r0 + k / ntc) as u32, (k % ntc) as u32);
                let (r0, c0) = (ti * tile_nrow, tj * tile_ncol);
                let (r1, c1) = ((r0 + tile_nrow).min(grid.nrow), (c0 + tile_ncol).min(grid.ncol));
                // Outer edges come from the grid itself so the last tiles end
                // exactly at xmax and ymin.
                let xmax = if c1 == grid.ncol { grid.xmax } else { grid.xmin + c1 as f64 * dx };
                let ymin = if r1 == grid.nrow { grid.ymin } else { grid.ymax - r1 as f64 * dy };
                let tgrid =
                    GridSpec::new(grid.xmin + c0 as f64 * dx, ymin, xmax, grid.ymax - r0 as f64 * dy, c1 - c0, r1 - r0);
                Tile { tile_row: ti, tile_col: tj, row_off: r0, col_off: c0, grid: tgrid, result }
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{burn, BurnOptions, Coord, Geometry, Polygon};
    use proptest::prelude::*;

    fn ring(cx: f64, cy: f64, r: f64, k: usize, wobble: &[f64]) -> Vec<Coord> {
        let mut v: Vec<Coord> = (0..k)
            .map(|i| {
                let a = i as f64 / k as f64 * std::f64::consts::TAU;
                let rr = r * (0.6 + 0.4 * wobble[i % wobble.len()]);
                Coord::new(cx + rr * a.cos(), cy + rr * a.sin())
            })
            .collect();
        v.push(v[0]);
        v
    }

    proptest! {
        /// The parallel slice equals the serial one, tile for tile and
        /// record for record, whatever the threads, tile shape and input
        /// (the parallel path forced on from one record).
        #[test]
        fn parallel_equals_serial(
            centres in prop::collection::vec((0.0..37.0f64, 0.0..23.0f64, 0.3..6.0f64), 0..24),
            wobble in prop::collection::vec(0.0..1.0f64, 9),
            tc in 1u32..20, tr in 1u32..12, threads in 2usize..12, approx in any::<bool>(),
        ) {
            let grid = GridSpec::new(0.0, 0.0, 37.0, 23.0, 53, 31);
            let mut geoms: Vec<Geometry> = centres
                .iter()
                .map(|&(x, y, r)| Geometry::Polygon(Polygon::new(vec![ring(x, y, r, 9, &wobble)])))
                .collect();
            if let Some(&(x, y, r)) = centres.first() {
                geoms.push(Geometry::LineString(ring(y, x * 0.6, r, 7, &wobble)));
                geoms.push(Geometry::MultiPoint((0..53).map(|c| Coord::new(c as f64 * 37.0 / 53.0, (c % 23) as f64)).collect()));
            }
            let opts = if approx { BurnOptions::approx() } else { BurnOptions::default() };
            let r = burn(&geoms, &grid, opts).unwrap();
            let want = slice_tiles(&r, &grid, tc, tr);
            let got = slice_tiles_par_from(&r, &grid, tc, tr, threads, 0);
            prop_assert_eq!(got, want);
        }
    }

    #[test]
    #[should_panic(expected = "outside the 4 x 4 grid")]
    fn parallel_panics_as_serial_on_a_record_outside_the_grid() {
        let grid = GridSpec::new(0.0, 0.0, 4.0, 4.0, 4, 4);
        let mut r = BurnResult::default();
        r.edges.push(GridEdge { row: 1, col: 1, fraction: 0.5, id: 0 });
        r.edges.push(GridEdge { row: 9, col: 1, fraction: 0.5, id: 0 });
        slice_tiles_par_from(&r, &grid, 2, 2, 4, 0);
    }

    #[test]
    fn falls_back_to_serial() {
        // one thread, one tile row, or few records: same answer, no threads
        let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
        let sq = Geometry::Polygon(Polygon::new(vec![ring(5.0, 5.0, 3.0, 12, &[0.5])]));
        let r = burn(&[sq], &grid, BurnOptions::default()).unwrap();
        for (tc, tr, threads) in [(3, 3, 1), (3, 10, 8), (3, 3, 0), (3, 3, 8)] {
            assert_eq!(slice_tiles_par(&r, &grid, tc, tr, threads), slice_tiles(&r, &grid, tc, tr));
        }
    }
}
