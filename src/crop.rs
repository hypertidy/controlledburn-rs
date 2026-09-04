//! Sparse crop to a target window.
//!
//! [`crop`] is the Rust analogue of the R package's `crop_burn()` and the
//! Python `BurnResult.crop()`: it filters and clips the four sparse tables
//! (runs, edges, lines, points) to a sub-window of the grid, re-basing all
//! row/col indices to 1 and snapping the window outward to whole cells. No
//! dense buffer is allocated — this is structured filtering only, and it
//! chains straight into [`crate::materialize`] to cut a tile from a single
//! burn.
//!
//! Unlike the R/Python bindings, a [`BurnResult`] carries no grid metadata,
//! so the grid the result was burned into is passed explicitly and the
//! snapped sub-grid is returned alongside the cropped tables.

use crate::grid::GridSpec;
use crate::output::BurnResult;

/// Crop `result` (burned into `grid`) to a target window, returning the
/// cropped tables together with the snapped sub-grid.
///
/// `target` is `[xmin, ymin, xmax, ymax]` in CRS units, matching the field
/// order of [`GridSpec`]. (Note this differs from the R/Python bindings,
/// whose `target` is `[xmin, xmax, ymin, ymax]`.) The window is snapped
/// outward to cell boundaries and clamped to the grid.
///
/// Returns `None` if the target does not overlap the grid — the idiomatic
/// counterpart to the bindings' warn-and-return-empty behaviour, which
/// would otherwise require a zero-dimension grid that [`GridSpec::validate`]
/// rejects. `notes` are carried through unchanged, since they are indexed by
/// input geometry rather than by cell.
///
/// # Examples
///
/// ```
/// use controlledburn::{burn, crop, materialize, BurnOptions, Coord,
///     Geometry, GridSpec, MaterializeOptions, Polygon};
///
/// let square = Geometry::Polygon(Polygon::new(vec![vec![
///     Coord::new(2.5, 4.5), Coord::new(6.5, 4.5), Coord::new(6.5, 8.5),
///     Coord::new(2.5, 8.5), Coord::new(2.5, 4.5),
/// ]]));
/// let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
/// let r = burn(&[square], &grid, BurnOptions::default()).unwrap();
///
/// let (sub, subgrid) = crop(&r, &grid, [3.0, 5.0, 7.0, 9.0]).unwrap();
/// assert_eq!((subgrid.ncol, subgrid.nrow), (4, 4));
///
/// let mut tile = vec![f64::NAN; (subgrid.ncol * subgrid.nrow) as usize];
/// materialize(&sub, &mut tile, subgrid.ncol, subgrid.nrow, None,
///     &MaterializeOptions::default()).unwrap();
/// ```
pub fn crop(result: &BurnResult, grid: &GridSpec, target: [f64; 4]) -> Option<(BurnResult, GridSpec)> {
    let [txmin, tymin, txmax, tymax] = target;
    let dx = grid.dx();
    let dy = grid.dy();
    let ncol = grid.ncol as i64;
    let nrow = grid.nrow as i64;

    // 1-based inclusive limits, snapped outward. The eps nudge avoids
    // floor/ceil flips when the target aligns exactly with a cell boundary.
    const EPS: f64 = 1e-8;
    let col_lo = ((((txmin - grid.xmin) / dx) + EPS).floor() as i64 + 1).max(1);
    let col_hi = ((((txmax - grid.xmin) / dx) - EPS).ceil() as i64).min(ncol);
    let row_hi = ((((grid.ymax - tymin) / dy) - EPS).ceil() as i64).min(nrow);
    let row_lo = ((((grid.ymax - tymax) / dy) + EPS).floor() as i64 + 1).max(1);

    if col_lo > col_hi || row_lo > row_hi {
        return None;
    }

    let (col_lo, col_hi) = (col_lo as i32, col_hi as i32);
    let (row_lo, row_hi) = (row_lo as i32, row_hi as i32);

    let subgrid = GridSpec::new(
        grid.xmin + (col_lo - 1) as f64 * dx,
        grid.ymax - row_hi as f64 * dy,
        grid.xmin + col_hi as f64 * dx,
        grid.ymax - (row_lo - 1) as f64 * dy,
        (col_hi - col_lo + 1) as u32,
        (row_hi - row_lo + 1) as u32,
    );

    let in_rows = |row: i32| row >= row_lo && row <= row_hi;
    let in_cols = |col: i32| col >= col_lo && col <= col_hi;

    let mut out = BurnResult::default();

    for r in &result.runs {
        // A run overlaps the window when its span meets [col_lo, col_hi].
        if in_rows(r.row) && r.col_end >= col_lo && r.col_start <= col_hi {
            out.runs.push(crate::output::GridRun {
                row: r.row - row_lo + 1,
                col_start: r.col_start.max(col_lo) - col_lo + 1,
                col_end: r.col_end.min(col_hi) - col_lo + 1,
                id: r.id,
            });
        }
    }
    for e in &result.edges {
        if in_rows(e.row) && in_cols(e.col) {
            out.edges.push(crate::output::GridEdge {
                row: e.row - row_lo + 1,
                col: e.col - col_lo + 1,
                fraction: e.fraction,
                id: e.id,
            });
        }
    }
    for l in &result.lines {
        if in_rows(l.row) && in_cols(l.col) {
            out.lines.push(crate::output::GridLine {
                row: l.row - row_lo + 1,
                col: l.col - col_lo + 1,
                length: l.length,
                id: l.id,
            });
        }
    }
    for p in &result.points {
        if in_rows(p.row) && in_cols(p.col) {
            out.points.push(crate::output::GridPoint { row: p.row - row_lo + 1, col: p.col - col_lo + 1, id: p.id });
        }
    }
    out.notes = result.notes.clone();

    Some((out, subgrid))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::burn::{burn, BurnOptions};
    use crate::geometry::{Coord, Geometry, Polygon};
    use crate::materialize::{materialize, EdgePolicy, MaterializeOptions, PixelFn};

    // extent [0,10] x [0,10], 10x10 cells (cell size 1).
    fn grid10() -> GridSpec {
        GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10)
    }

    // shapely.box(2.5, 4.5, 6.5, 8.5)
    fn box_poly(xmin: f64, ymin: f64, xmax: f64, ymax: f64) -> Geometry {
        Geometry::Polygon(Polygon::new(vec![vec![
            Coord::new(xmin, ymin),
            Coord::new(xmax, ymin),
            Coord::new(xmax, ymax),
            Coord::new(xmin, ymax),
            Coord::new(xmin, ymin),
        ]]))
    }

    fn box_result() -> BurnResult {
        burn(&[box_poly(2.5, 4.5, 6.5, 8.5)], &grid10(), BurnOptions::default()).unwrap()
    }

    // materialize with sum + fraction over a NaN background; returns a dense
    // row-major buffer, treating never-touched cells as NaN.
    fn dense(r: &BurnResult, g: &GridSpec) -> Vec<f64> {
        let mut buf = vec![f64::NAN; (g.ncol * g.nrow) as usize];
        let opts = MaterializeOptions { fn_: PixelFn::Sum, edge_policy: EdgePolicy::Fraction, threshold: 0.5 };
        materialize(r, &mut buf, g.ncol, g.nrow, None, &opts).unwrap();
        buf
    }

    // NaN-aware near-equality for dense buffers (NaN == NaN as background).
    fn approx_eq(a: &[f64], b: &[f64]) {
        assert_eq!(a.len(), b.len(), "buffer lengths differ");
        for (i, (x, y)) in a.iter().zip(b).enumerate() {
            match (x.is_nan(), y.is_nan()) {
                (true, true) => {}
                (false, false) => assert!((x - y).abs() < 1e-9, "cell {i}: {x} vs {y}"),
                _ => panic!("cell {i}: NaN mismatch {x} vs {y}"),
            }
        }
    }

    // Slice rows [r0,r1] x cols [c0,c1] (1-based inclusive) from a full
    // row-major buffer with `ncol` columns.
    fn slice(full: &[f64], ncol: i32, r0: i32, r1: i32, c0: i32, c1: i32) -> Vec<f64> {
        let mut out = Vec::new();
        for r in r0..=r1 {
            for c in c0..=c1 {
                out.push(full[((r - 1) * ncol + (c - 1)) as usize]);
            }
        }
        out
    }

    #[test]
    fn overlap_shape_and_extent() {
        let (_, g) = crop(&box_result(), &grid10(), [3.0, 5.0, 7.0, 9.0]).unwrap();
        assert_eq!((g.nrow, g.ncol), (4, 4));
        assert_eq!((g.xmin, g.ymin, g.xmax, g.ymax), (3.0, 5.0, 7.0, 9.0));
    }

    #[test]
    fn rebases_indices_to_one() {
        let (sub, g) = crop(&box_result(), &grid10(), [3.0, 5.0, 7.0, 9.0]).unwrap();
        for r in &sub.runs {
            assert!(r.row >= 1 && r.row <= g.nrow as i32);
            assert!(r.col_start >= 1 && r.col_end <= g.ncol as i32);
        }
        for e in &sub.edges {
            assert!(e.row >= 1 && e.row <= g.nrow as i32);
            assert!(e.col >= 1 && e.col <= g.ncol as i32);
        }
    }

    #[test]
    fn crop_matches_full_slice() {
        let r = box_result();
        let full = dense(&r, &grid10());
        let (sub, g) = crop(&r, &grid10(), [3.0, 5.0, 7.0, 9.0]).unwrap();
        // window x in [3,7], y in [5,9] -> rows 2..5, cols 4..7 (1-based)
        approx_eq(&dense(&sub, &g), &slice(&full, 10, 2, 5, 4, 7));
    }

    #[test]
    fn edge_aligned_window() {
        let r = box_result();
        let full = dense(&r, &grid10());
        let (sub, g) = crop(&r, &grid10(), [2.0, 4.0, 6.0, 8.0]).unwrap();
        assert_eq!((g.nrow, g.ncol), (4, 4));
        assert_eq!((g.xmin, g.ymin, g.xmax, g.ymax), (2.0, 4.0, 6.0, 8.0));
        // y in [4,8] -> rows 3..6, x in [2,6] -> cols 3..6 (1-based)
        approx_eq(&dense(&sub, &g), &slice(&full, 10, 3, 6, 3, 6));
    }

    #[test]
    fn full_window_is_identity() {
        let r = box_result();
        let (sub, g) = crop(&r, &grid10(), [0.0, 0.0, 10.0, 10.0]).unwrap();
        assert_eq!((g.nrow, g.ncol), (10, 10));
        assert_eq!((g.xmin, g.ymin, g.xmax, g.ymax), (0.0, 0.0, 10.0, 10.0));
        approx_eq(&dense(&sub, &g), &dense(&r, &grid10()));
    }

    #[test]
    fn window_larger_than_grid_clamps() {
        let r = box_result();
        let (sub, g) = crop(&r, &grid10(), [-5.0, -5.0, 15.0, 15.0]).unwrap();
        assert_eq!((g.nrow, g.ncol), (10, 10));
        assert_eq!((g.xmin, g.ymin, g.xmax, g.ymax), (0.0, 0.0, 10.0, 10.0));
        approx_eq(&dense(&sub, &g), &dense(&r, &grid10()));
    }

    #[test]
    fn clips_and_rebases_runs() {
        // polygon covering the whole grid: every row is one run, cols 1..10
        let r = burn(&[box_poly(0.0, 0.0, 10.0, 10.0)], &grid10(), BurnOptions::default()).unwrap();
        assert!(r.edges.is_empty()); // grid-aligned -> pure interior
        let full = dense(&r, &grid10());
        let (sub, g) = crop(&r, &grid10(), [3.0, 3.0, 7.0, 7.0]).unwrap();
        assert_eq!((g.nrow, g.ncol), (4, 4));
        assert!(sub.runs.iter().all(|x| x.col_start >= 1 && x.col_end <= 4));
        assert!(sub.runs.iter().all(|x| x.row >= 1 && x.row <= 4));
        approx_eq(&dense(&sub, &g), &slice(&full, 10, 4, 7, 4, 7));
    }

    #[test]
    fn non_overlap_returns_none() {
        assert!(crop(&box_result(), &grid10(), [100.0, 100.0, 200.0, 200.0]).is_none());
    }

    #[test]
    fn line_table_cropped_and_rebased() {
        let line = Geometry::LineString(vec![Coord::new(0.5, 0.5), Coord::new(9.5, 9.5)]);
        let r = burn(&[line], &grid10(), BurnOptions::default()).unwrap();
        assert!(!r.lines.is_empty());
        let (sub, g) = crop(&r, &grid10(), [3.0, 3.0, 7.0, 7.0]).unwrap();
        assert_eq!((g.nrow, g.ncol), (4, 4));
        assert!(!sub.lines.is_empty());
        assert!(sub.lines.len() <= r.lines.len());
        assert!(sub.lines.iter().all(|l| l.row >= 1 && l.row <= 4 && l.col >= 1 && l.col <= 4));
    }

    #[test]
    fn point_table_cropped_and_rebased() {
        let pts = vec![Geometry::Point(Coord::new(1.5, 1.5)), Geometry::Point(Coord::new(5.5, 5.5))];
        let r = burn(&pts, &grid10(), BurnOptions::default()).unwrap();
        let (sub, _) = crop(&r, &grid10(), [4.0, 4.0, 8.0, 8.0]).unwrap();
        // only (5.5, 5.5) falls inside the window
        assert_eq!(sub.points.len(), 1);
        assert!(sub.points.iter().all(|p| p.row >= 1 && p.row <= 4 && p.col >= 1 && p.col <= 4));
    }
}
