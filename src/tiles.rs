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

use crate::grid::GridSpec;
use crate::output::{BurnResult, GridEdge, GridLine, GridPoint, GridRun};

/// One non-empty tile of a [`slice_tiles`] call.
#[derive(Clone, Debug)]
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
    assert!(tile_ncol > 0 && tile_nrow > 0, "tile shape must be at least 1 x 1");
    let (tc, tr) = (tile_ncol as i32, tile_nrow as i32);
    let ntc = (grid.ncol as usize + tile_ncol as usize - 1) / tile_ncol as usize;
    let ntr = (grid.nrow as usize + tile_nrow as usize - 1) / tile_nrow as usize;
    let at = |row: i32, col: i32| -> usize {
        assert!(
            row >= 0 && col >= 0 && (row as u32) < grid.nrow && (col as u32) < grid.ncol,
            "record ({row}, {col}) is outside the {} x {} grid",
            grid.ncol,
            grid.nrow
        );
        (row / tr) as usize * ntc + (col / tc) as usize
    };

    let mut out: Vec<BurnResult> = vec![BurnResult::default(); ntc * ntr];
    for run in &result.runs {
        let lr = run.row % tr;
        let mut c = run.col_start;
        while c < run.col_end {
            let c0 = (c / tc) * tc;
            let ce = (c0 + tc).min(run.col_end);
            out[at(run.row, c)].runs.push(GridRun { row: lr, col_start: c - c0, col_end: ce - c0, id: run.id });
            c = ce;
        }
    }
    for e in &result.edges {
        out[at(e.row, e.col)].edges.push(GridEdge { row: e.row % tr, col: e.col % tc, ..*e });
    }
    for l in &result.lines {
        out[at(l.row, l.col)].lines.push(GridLine { row: l.row % tr, col: l.col % tc, ..*l });
    }
    for p in &result.points {
        out[at(p.row, p.col)].points.push(GridPoint { row: p.row % tr, col: p.col % tc, ..*p });
    }

    let (dx, dy) = (grid.dx(), grid.dy());
    out.into_iter()
        .enumerate()
        .filter(|(_, r)| !(r.runs.is_empty() && r.edges.is_empty() && r.lines.is_empty() && r.points.is_empty()))
        .map(|(k, result)| {
            let (ti, tj) = ((k / ntc) as u32, (k % ntc) as u32);
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
