//! The cell-stepping walker: follows a polyline through the padded grid,
//! recording per-cell traversals (entry side, coordinates, exit side).
//!
//! This is a port of controlledburn's `walk_polyline`, itself a rewrite of
//! exactextract's `RasterCellIntersection::process_line` without the
//! `Cell` class. The control flow is kept identical, including the
//! "incomplete initial traversal" re-append for closed rings, so that the
//! set of traversals per cell (and hence every coverage fraction) matches
//! the C++ bit for bit.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

use crate::ee::bbox::BBox;
use crate::ee::side::Side;
use crate::geometry::Coord;
use crate::grid::Grid;

/// One pass of the polyline through one cell. Coordinates live in the
/// walk's shared arena (`Cells::coords`) as the range `start..start + len`.
///
/// 16 bytes: a ring of 10^6 vertices can make 10^8 traversals.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Traversal {
    start: u32,
    len: u32,
    /// the next traversal of the same cell in `Cells::flat`, or `NONE`
    next: u32,
    pub entry: Option<Side>,
    pub exit: Option<Side>,
}

const NONE: u32 = u32::MAX;

impl Traversal {
    #[inline]
    pub fn coords<'a>(&self, arena: &'a [Coord]) -> &'a [Coord] {
        let s = self.start as usize;
        &arena[s..s + self.len as usize]
    }

    #[inline]
    pub fn traversed(&self) -> bool {
        self.entry.is_some() && self.exit.is_some()
    }

    #[inline]
    pub fn is_closed_ring(&self, arena: &[Coord]) -> bool {
        let c = self.coords(arena);
        c.len() >= 3 && c[0] == c[c.len() - 1]
    }

    pub fn multiple_unique_coordinates(&self, arena: &[Coord]) -> bool {
        let c = self.coords(arena);
        c[1..].iter().any(|x| *x != c[0])
    }
}

/// Per-cell traversal data for one cell (row, col) of the padded grid.
/// Its traversals are `Cells::traversals(self)`.
#[derive(Clone, Debug)]
pub(crate) struct CellRecord {
    pub row: usize,
    pub col: usize,
    pub bbox: BBox,
    /// first and last of its traversals in `Cells::flat`
    head: u32,
    tail: u32,
}

/// Cells visited by a walk, in first-visit order, plus the coordinate
/// arena all traversals index into. Consumers that need a deterministic
/// order sort by (row, col).
///
/// Every traversal lives in one array, in walk order, each linked to the
/// next traversal of its cell, rather than in a vector per cell: a
/// detailed ring of 10^6 vertices visits ten million cells about eleven
/// times each, and the vectors were 27 million allocations.
/// Grouping the array by cell after the walk would make the reads
/// contiguous, but costs more than it saves (an in-place counting sort
/// took 7 s on that ring) or a second copy of the array.
#[derive(Clone, Debug, Default)]
pub(crate) struct Cells {
    pub records: Vec<CellRecord>,
    pub coords: Vec<Coord>,
    flat: Vec<Traversal>,
}

impl Cells {
    /// The traversals of `cr`, in the order the walk made them.
    #[inline]
    pub fn traversals<'a>(&'a self, cr: &CellRecord) -> Traversals<'a> {
        Traversals { flat: &self.flat, at: cr.head }
    }
}

/// The traversals of one cell, in walk order (`Cells::traversals`).
pub(crate) struct Traversals<'a> {
    flat: &'a [Traversal],
    at: u32,
}

impl<'a> Iterator for Traversals<'a> {
    type Item = &'a Traversal;

    #[inline]
    fn next(&mut self) -> Option<&'a Traversal> {
        if self.at == NONE {
            return None;
        }
        let t = &self.flat[self.at as usize];
        self.at = t.next;
        Some(t)
    }
}

/// Minimal multiplicative hasher for the (row, col) index. The walker
/// mostly revisits the current cell or a neighbour, so the map stays hot.
#[derive(Default)]
struct CellHasher(u64);

impl Hasher for CellHasher {
    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        for b in bytes {
            self.write_u64(*b as u64);
        }
    }
    #[inline]
    fn write_u64(&mut self, v: u64) {
        self.0 = (self.0.rotate_left(5) ^ v).wrapping_mul(0x517c_c1b7_2722_0a95);
    }
    #[inline]
    fn write_usize(&mut self, v: usize) {
        self.write_u64(v as u64);
    }
}

type CellIndex = HashMap<(usize, usize), usize, BuildHasherDefault<CellHasher>>;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Location {
    Inside,
    Boundary,
    Outside,
}

#[inline]
fn point_location(bbox: &BBox, c: &Coord) -> Location {
    if bbox.strictly_contains(c) {
        Location::Inside
    } else if bbox.contains(c) {
        Location::Boundary
    } else {
        Location::Outside
    }
}

/// Walk `coords` through `grid`, returning the visited cells.
///
/// `closed` selects how an initial traversal that began strictly inside a
/// cell (no entry side) is handled once it exits:
///   - `true` (polygon ring): its coordinates are appended to the input so
///     the start cell is re-entered with a proper entry side on the second
///     pass, stitching the start back to the closing edge.
///   - `false` (open polyline): the partial traversal is kept as is; the
///     line consumer accepts `entry == None` as valid.
///
/// Errors correspond to the C++ `std::runtime_error` throws ("Never get
/// here") and are reported by the caller as a note.
/// Snap vertices that lie within a tiny tolerance of a cell edge exactly
/// onto that edge, using the same arithmetic `Grid::cell` uses to place
/// the edges.
///
/// A vertex parsed from input (0.3) and an edge computed from the grid
/// (3 * 0.1 = 0.30000000000000004) can differ by a few ulp. Left unsnapped,
/// such a vertex is classified OUTSIDE the cell it visually sits on the
/// corner of, `BBox::crossing` extrapolates past the segment end, and the
/// walk leaves the grid along the extended line (the C++ core runs until
/// `bad_alloc`; here the padded-grid bounds check turns it into an error
/// and an empty polygon). Where the walk does close, a degenerate corner
/// exit makes `left_hand_area` over-count the cell. Tolerance is relative
/// to the cell size; only degenerate near-edge vertices move, by less than
/// 1e-9 of a cell.
fn snap_to_lattice(coords: &mut [Coord], grid: &Grid<1>) {
    let (dx, dy) = (grid.dx(), grid.dy());
    if dx <= 0.0 || dy <= 0.0 || dx.is_nan() || dy.is_nan() {
        return;
    }
    let (x0, x1, y0, y1) = (grid.xmin(), grid.xmax(), grid.ymin(), grid.ymax());
    // padded grids carry one padding row/col on each side
    let ncol = (grid.cols() - 2) as f64;
    let nrow = (grid.rows() - 2) as f64;
    let (tol_x, tol_y) = (dx * 1e-9, dy * 1e-9);
    for c in coords.iter_mut() {
        let kx = ((c.x - x0) / dx).round();
        if kx >= 0.0 && kx <= ncol {
            let bx = if kx == ncol { x1 } else { x0 + kx * dx };
            if (c.x - bx).abs() <= tol_x {
                c.x = bx;
            }
        }
        let ky = ((y1 - c.y) / dy).round();
        if ky >= 0.0 && ky <= nrow {
            let by = if ky == nrow { y0 } else { y1 - ky * dy };
            if (c.y - by).abs() <= tol_y {
                c.y = by;
            }
        }
    }
}

/// Move `c` (on the boundary of `bbox`) onto the corner of `bbox` it lies
/// within `tol_x`, `tol_y` of, if any. Corners of padding cells are at
/// `f64::MIN`/`f64::MAX` and are never within tolerance.
#[inline]
fn snap_to_corner(c: Coord, bbox: &BBox, tol_x: f64, tol_y: f64) -> Coord {
    let x = if (c.x - bbox.xmin).abs() <= tol_x {
        bbox.xmin
    } else if (c.x - bbox.xmax).abs() <= tol_x {
        bbox.xmax
    } else {
        return c;
    };
    let y = if (c.y - bbox.ymin).abs() <= tol_y {
        bbox.ymin
    } else if (c.y - bbox.ymax).abs() <= tol_y {
        bbox.ymax
    } else {
        return c;
    };
    Coord::new(x, y)
}

pub(crate) fn walk_polyline(coords: Vec<Coord>, grid: &Grid<1>, closed: bool) -> Result<Cells, String> {
    let mut w = Walker { input: coords, ..Walker::default() };
    w.walk(grid, closed)?;
    Ok(w.cells)
}

/// The walk's buffers, kept between walks so that burning many small
/// rings does not allocate per ring and per cell: the input copy, the
/// visited cells and their arena, the cell index, and the per-cell
/// traversal lists (cleared and handed back to the next walk's cells).
#[derive(Default)]
pub(crate) struct Walker {
    /// The coordinates to walk. The caller fills it; the walk snaps it
    /// and may extend it (the re-append for closed rings).
    pub input: Vec<Coord>,
    /// The cells of the last walk.
    pub cells: Cells,
    index: CellIndex,
}

impl Walker {
    fn reset(&mut self) {
        self.cells.records.clear();
        self.cells.flat.clear();
        self.cells.coords.clear();
        // Clearing a hash map costs its capacity, not its length: after one
        // large ring, a run of small ones would each pay for it. Start
        // afresh once the map is much larger than what the last walk used.
        if self.index.capacity() > 4 * self.index.len().max(64) {
            self.index = CellIndex::default();
        } else {
            self.index.clear();
        }
    }

    /// Walk `self.input` through `grid`; the visited cells are left in
    /// `self.cells`. See `walk_polyline`.
    pub fn walk(&mut self, grid: &Grid<1>, closed: bool) -> Result<(), String> {
        self.reset();
        let Walker { input: coords, cells, index } = self;
        if coords.is_empty() {
            return Ok(());
        }
        snap_to_lattice(coords, grid);
        let n_input = coords.len();
        // Every input vertex produces at most one arena entry plus one exit
        // crossing per cell change; reserve generously once.
        cells.coords.reserve(coords.len() * 2 + 8);
        let arena = &mut cells.coords;
        let records = &mut cells.records;
        let flat = &mut cells.flat;
        // Cache of the last cell looked up: consecutive visits are usually
        // the same cell or a neighbour.
        let mut last: Option<((usize, usize), usize)> = None;
        let n_rows = grid.rows();
        let n_cols = grid.cols();
        // Same tolerance as `snap_to_lattice`, from the grid's cell size (a
        // padding cell's box is unbounded).
        let (tol_x, tol_y) = (grid.dx() * 1e-9, grid.dy() * 1e-9);

        let mut pos: usize = 0;
        let mut row = grid.get_row(coords[0].y);
        let mut col = grid.get_column(coords[0].x);
        let mut last_exit: Option<Coord> = None;

        while pos < coords.len() {
            if row >= n_rows || col >= n_cols {
                return Err("walker stepped outside the padded grid".into());
            }
            let ci = match last {
                Some((k, i)) if k == (row, col) => i,
                _ => {
                    let i = *index.entry((row, col)).or_insert_with(|| {
                        records.push(CellRecord { row, col, bbox: grid.cell(row, col), head: NONE, tail: NONE });
                        records.len() - 1
                    });
                    last = Some(((row, col), i));
                    i
                }
            };
            let bbox = records[ci].bbox;

            // 32-bit offsets: an arena of 2^32 coordinates is 64 GB
            let start = u32::try_from(arena.len()).map_err(|_| "ring walk too large (2^32 coordinates)")?;
            let mut trav = Traversal { start, len: 0, next: NONE, entry: None, exit: None };

            while pos < coords.len() {
                let next = match last_exit {
                    Some(c) => c,
                    None => coords[pos],
                };

                if trav.len == 0 {
                    // First coordinate for this traversal: enter the cell.
                    trav.entry = bbox.side(&next);
                    arena.push(next);
                    trav.len += 1;
                    if last_exit.is_some() {
                        last_exit = None;
                    } else {
                        pos += 1;
                    }
                    continue;
                }

                if point_location(&bbox, &next) != Location::Outside {
                    arena.push(next);
                    trav.len += 1;
                    if last_exit.is_some() {
                        last_exit = None;
                    } else {
                        pos += 1;
                    }
                } else {
                    // Outside: compute the exit crossing. Use the previous
                    // ORIGINAL vertex for robustness (same as Cell::take).
                    let from = if pos > 0 { coords[pos - 1] } else { arena[arena.len() - 1] };
                    let x = bbox.crossing(&from, &next).ok_or_else(|| String::from("Never get here."))?;
                    let (mut exit_c, mut exit_s) = (x.coord, x.side);
                    // Rounding guard (see snap_to_lattice): a crossing that lies
                    // beyond the segment end can only come from a target that is
                    // outside the box by rounding error. Clamp the target onto
                    // the box and exit there instead of extrapolating.
                    {
                        let d_next = (next.x - from.x).powi(2) + (next.y - from.y).powi(2);
                        let d_exit = (exit_c.x - from.x).powi(2) + (exit_c.y - from.y).powi(2);
                        if !exit_c.x.is_finite() || !exit_c.y.is_finite() || d_exit > d_next * (1.0 + 1e-9) {
                            exit_c = Coord::new(next.x.clamp(bbox.xmin, bbox.xmax), next.y.clamp(bbox.ymin, bbox.ymax));
                            exit_s = if next.x < bbox.xmin {
                                Side::Left
                            } else if next.x > bbox.xmax {
                                Side::Right
                            } else if next.y < bbox.ymin {
                                Side::Bottom
                            } else {
                                Side::Top
                            };
                        }
                    }
                    // Corner guard: an edge through a cell corner gives a
                    // crossing an ulp or two off the corner, on the wall
                    // either side of it. The next cell then starts a sliver
                    // traversal from there to the corner, or from the corner
                    // to a point an ulp down the wall, where the true path
                    // only touches the corner. Its two ends can be out of
                    // order around the cell boundary, which the old chain
                    // chaser read as a slab of the cell or all of it
                    // (controlledburn-rs#6). Exit through the corner itself;
                    // the sliver becomes a single-point touch, which the
                    // coverage pass already ignores.
                    exit_c = snap_to_corner(exit_c, &bbox, tol_x, tol_y);
                    arena.push(exit_c);
                    trav.len += 1;
                    trav.exit = Some(exit_s);
                    if exit_c != next {
                        last_exit = Some(exit_c);
                    }
                    break;
                }
            }

            // Force exit if stuck on the boundary (Cell::force_exit).
            if trav.exit.is_none() && trav.len > 0 {
                let last = arena[arena.len() - 1];
                if point_location(&bbox, &last) == Location::Boundary {
                    trav.exit = bbox.side(&last);
                }
            }

            let exited = trav.exit.is_some();
            let incomplete = exited && trav.entry.is_none();

            if incomplete && closed {
                // A closed ring is reprocessed at most once per starting cell;
                // if the coordinate list keeps growing the walk is not converging.
                if coords.len() > 4 * n_input + 64 {
                    return Err("ring walk did not converge (vertex on cell boundary?)".into());
                }
                coords.extend_from_slice(trav.coords(arena));
            }

            let exit_side = trav.exit;
            let k =
                u32::try_from(flat.len()).ok().filter(|&k| k != NONE).ok_or("ring walk too large (2^32 traversals)")?;
            let rec = &mut records[ci];
            if rec.head == NONE {
                rec.head = k;
            } else {
                flat[rec.tail as usize].next = k;
            }
            rec.tail = k;
            flat.push(trav);

            if exited {
                match exit_side {
                    Some(Side::Top) => row = row.checked_sub(1).ok_or("walker stepped above the padded grid")?,
                    Some(Side::Bottom) => row += 1,
                    Some(Side::Left) => col = col.checked_sub(1).ok_or("walker stepped left of the padded grid")?,
                    Some(Side::Right) => col += 1,
                    None => {}
                }
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::GridSpec;

    #[test]
    fn diagonal_line_touches_expected_cells() {
        let g = GridSpec::new(0., 0., 4., 4., 4, 4).bounded().make_infinite();
        let line = vec![Coord::new(0.5, 0.4), Coord::new(3.5, 3.4)];
        let cells = walk_polyline(line, &g, false).unwrap();
        // A 45-degree line offset from the corners clips 7 cells
        // (padded indices: row 4..1, col 1..4).
        let mut keys: Vec<_> = cells.records.iter().map(|c| (c.row, c.col)).collect();
        keys.sort();
        assert_eq!(keys, vec![(1, 4), (2, 3), (2, 4), (3, 2), (3, 3), (4, 1), (4, 2)]);
        // start cell has entry None
        assert!(cells.traversals(&cells.records[0]).next().unwrap().entry.is_none());
    }

    #[test]
    fn closed_ring_reappends_start() {
        let g = GridSpec::new(0., 0., 10., 10., 10, 10).bounded().make_infinite();
        let ring = vec![
            Coord::new(2.5, 2.5),
            Coord::new(6.5, 2.5),
            Coord::new(6.5, 6.5),
            Coord::new(2.5, 6.5),
            Coord::new(2.5, 2.5),
        ];
        let cells = walk_polyline(ring, &g, true).unwrap();
        // The start cell (row 8, col 3 in padded coords) gets two
        // traversals: the incomplete one and the stitched one.
        let start = cells.records.iter().find(|c| (c.row, c.col) == (8, 3)).unwrap();
        let ts: Vec<_> = cells.traversals(start).collect();
        assert_eq!(ts.len(), 2);
        assert!(ts[0].entry.is_none());
        assert!(ts[1].traversed());
        assert_eq!(cells.records.len(), 16);
    }
}
