//! Cheap pre-pass capacity estimates for the output tables.
//!
//! Without reserving, every table grows by `Vec` doubling: a large burn
//! (hundreds of millions of runs) recopies gigabytes through reallocations
//! and faults each freshly committed page, which shows up as system time.
//! One extra pass over the input, reserving from these estimates, replaces
//! the doubling cascade with a single allocation of about the right size.
//!
//! The numbers are estimates, not hard bounds. Rows-spanned is a true
//! upper bound on a convex polygon's runs but only an estimate for a
//! non-convex one (several spans per row); the boundary-cell count is an
//! upper-ish estimate. Being off by a small factor still collapses ~28
//! doublings to one allocation, which is the whole point.

use crate::burn::BurnMode;
use crate::geometry::{Coord, Geometry};
use crate::grid::GridSpec;
use crate::output::Capacity;

/// Grid cells a polyline passes through, upper-ish: a segment crossing `v`
/// vertical and `h` horizontal grid lines visits `v + h + 1` cells.
fn polyline_cells(seq: &[Coord], dx: f64, dy: f64) -> usize {
    if seq.len() < 2 {
        return seq.len();
    }
    let mut n = 0.0f64;
    for w in seq.windows(2) {
        n += (w[1].x - w[0].x).abs() / dx + (w[1].y - w[0].y).abs() / dy;
    }
    if !n.is_finite() {
        return 0;
    }
    (n.ceil() as usize).saturating_add(1)
}

/// Number of grid rows a y-range spans, clamped to the grid.
fn rows_spanned(ymin: f64, ymax: f64, gs: &GridSpec) -> usize {
    let lo = ymin.max(gs.ymin);
    let hi = ymax.min(gs.ymax);
    if hi <= lo {
        return 0;
    }
    let rows = ((hi - lo) / gs.dy()).ceil() as usize;
    rows.clamp(1, (gs.nrow as usize).max(1))
}

/// Estimate output-table capacities for `geoms` on `grid` in `mode`.
///
/// One pass over the coordinates, no parsing and no allocation. Feed the
/// result to [`crate::BurnResult::with_capacity`]; [`crate::burn`] does this
/// automatically.
pub fn estimate_capacity(geoms: &[Geometry], grid: &GridSpec, mode: BurnMode) -> Capacity {
    estimate_capacity_iter(geoms, grid, mode)
}

/// Estimate over any iterator of geometries. `burn_wkb` uses this because
/// it buffers parsed geometries rather than holding a slice.
pub(crate) fn estimate_capacity_iter<'a>(
    geoms: impl IntoIterator<Item = &'a Geometry>,
    grid: &GridSpec,
    mode: BurnMode,
) -> Capacity {
    let (dx, dy) = (grid.dx(), grid.dy());
    let mut cap = Capacity::default();

    for g in geoms {
        // Non-finite geometries are skipped by the burn; skip them here too,
        // so a bogus coordinate cannot inflate the estimate to usize::MAX.
        if !g.is_finite() {
            continue;
        }
        for poly in g.polygons() {
            let Some(ext) = poly.exterior() else { continue };
            if ext.len() < 4 {
                continue;
            }
            let (mut xmin, mut xmax, mut ymin, mut ymax) = (ext[0].x, ext[0].x, ext[0].y, ext[0].y);
            for c in ext {
                xmin = xmin.min(c.x);
                xmax = xmax.max(c.x);
                ymin = ymin.min(c.y);
                ymax = ymax.max(c.y);
            }
            // Skip polygons whose envelope misses the grid entirely.
            if xmax < grid.xmin || xmin > grid.xmax || ymax < grid.ymin || ymin > grid.ymax {
                continue;
            }
            // Interior runs: roughly one span per row the polygon spans.
            cap.runs += rows_spanned(ymin, ymax, grid);
            // Coverage boundary cells become edges; approx emits none.
            if mode == BurnMode::Coverage {
                for ring in &poly.rings {
                    cap.edges += polyline_cells(ring, dx, dy);
                }
            }
        }
        for line in g.lines() {
            cap.lines += polyline_cells(line, dx, dy);
        }
        cap.points += g.points().len();
    }
    cap
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::burn::{burn, BurnOptions};
    use crate::geometry::{Coord, Geometry, Polygon};
    use crate::output::BurnResult;

    fn square() -> Geometry {
        Geometry::Polygon(Polygon::new(vec![vec![
            Coord::new(2.5, 4.5),
            Coord::new(6.5, 4.5),
            Coord::new(6.5, 8.5),
            Coord::new(2.5, 8.5),
            Coord::new(2.5, 4.5),
        ]]))
    }

    #[test]
    fn with_capacity_reserves() {
        let out = BurnResult::with_capacity(Capacity { runs: 100, edges: 50, lines: 10, points: 5 });
        assert!(out.runs.capacity() >= 100);
        assert!(out.edges.capacity() >= 50);
        assert!(out.lines.capacity() >= 10);
        assert!(out.points.capacity() >= 5);
        assert!(out.is_empty());
    }

    #[test]
    fn estimate_bounds_the_square() {
        let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
        let geoms = [square()];

        let cov = estimate_capacity(&geoms, &grid, BurnMode::Coverage);
        let r = burn(&geoms, &grid, BurnOptions::coverage()).unwrap();
        assert!(cov.runs >= r.runs.len(), "runs {} < actual {}", cov.runs, r.runs.len());
        assert!(cov.edges >= r.edges.len(), "edges {} < actual {}", cov.edges, r.edges.len());

        // Approx reserves runs but never edges.
        let app = estimate_capacity(&geoms, &grid, BurnMode::Approx);
        assert_eq!(app.edges, 0);
        assert!(app.runs > 0);
    }

    #[test]
    fn estimate_counts_lines_and_points() {
        let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
        let geoms = [
            Geometry::LineString(vec![Coord::new(0.5, 0.5), Coord::new(9.5, 7.5)]),
            Geometry::MultiPoint(vec![Coord::new(1.5, 1.5), Coord::new(8.5, 8.5)]),
        ];
        let cap = estimate_capacity(&geoms, &grid, BurnMode::Coverage);
        assert_eq!(cap.points, 2);
        let r = burn(&geoms, &grid, BurnOptions::default()).unwrap();
        assert!(cap.lines >= r.lines.len(), "lines {} < actual {}", cap.lines, r.lines.len());
    }

    #[test]
    fn polygon_outside_grid_contributes_nothing() {
        let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
        let far = Geometry::Polygon(Polygon::new(vec![vec![
            Coord::new(100.0, 100.0),
            Coord::new(101.0, 100.0),
            Coord::new(101.0, 101.0),
            Coord::new(100.0, 100.0),
        ]]));
        let cap = estimate_capacity(&[far], &grid, BurnMode::Coverage);
        assert_eq!(cap, Capacity::default());
    }
}
