//! Approx-mode polygon rasterization: the cell-centre rule (fasterize
//! semantics) via a direct edge/row-centre intersection sweep. No walker,
//! no coverage math. O(edges * rows touched).
//!
//! A cell is inside iff the winding number at its centre is non-zero,
//! with a left-inclusive convention: an edge crossing exactly at the cell
//! centre x counts as left of centre.

use crate::ee::bbox::BBox;
use crate::geometry::{is_ccw, Polygon};
use crate::grid::GridSpec;
use crate::output::GridRun;
use crate::polygon::ring_envelope;

#[derive(Clone, Copy)]
struct Intercept {
    x: f64,
    delta: i32,
}

pub(crate) fn process_polygon_approx(poly: &Polygon, gs: &GridSpec, poly_id: i32, runs: &mut Vec<GridRun>) {
    let Some(exterior) = poly.rings.first() else { return };
    if exterior.len() < 4 {
        return;
    }
    let dx = gs.dx();
    let dy = gs.dy();
    let ncol = gs.ncol as i64;
    let nrow = gs.nrow as i64;

    let env = ring_envelope(exterior);
    let grid_box = BBox::new(gs.xmin, gs.ymin, gs.xmax, gs.ymax);
    if !env.intersects(&grid_box) {
        return;
    }

    // Row range the polygon can touch (0-based, row 0 at the top).
    let row_lo = (((gs.ymax - env.ymax) / dy).floor() as i64).max(0);
    let row_hi = (((gs.ymax - env.ymin) / dy).floor() as i64).min(nrow - 1);
    if row_lo > row_hi {
        return;
    }
    let n_sweep_rows = (row_hi - row_lo + 1) as usize;
    let mut row_intercepts: Vec<Vec<Intercept>> = vec![Vec::new(); n_sweep_rows];

    for (ri, ring) in poly.rings.iter().enumerate() {
        if ring.len() < 4 {
            continue;
        }
        let exterior = ri == 0;
        let mut wf = if is_ccw(ring) { 1 } else { -1 };
        if !exterior {
            wf = -wf;
        }

        for w in ring.windows(2) {
            let (x0, y0) = (w[0].x, w[0].y);
            let (x1, y1) = (w[1].x, w[1].y);
            if y0 == y1 {
                continue; // horizontal edge
            }
            let (ya, yb) = if y0 > y1 { (y1, y0) } else { (y0, y1) };

            // Rows whose y_mid lies in (ya, yb]: top-inclusive,
            // bottom-exclusive, matching the walker's half-open crossing
            // convention. y_mid(r) = ymax - (r + 0.5) * dy.
            let e_row_lo = (((gs.ymax - yb) / dy - 0.5 - 1e-10).ceil() as i64).max(row_lo);
            let e_row_hi = ((((gs.ymax - ya) / dy - 0.5 - 1e-10).ceil() as i64) - 1).min(row_hi);
            if e_row_lo > e_row_hi {
                continue;
            }

            let inv_dy_edge = 1.0 / (y1 - y0);
            for r in e_row_lo..=e_row_hi {
                let y_mid = gs.ymax - (r as f64 + 0.5) * dy;
                let t = (y_mid - y0) * inv_dy_edge;
                let x_int = x0 + t * (x1 - x0);
                let delta = if y0 >= y_mid { -1 } else { 1 } * wf;
                row_intercepts[(r - row_lo) as usize].push(Intercept { x: x_int, delta });
            }
        }
    }

    for r in row_lo..=row_hi {
        let intercepts = &mut row_intercepts[(r - row_lo) as usize];
        if intercepts.is_empty() {
            continue;
        }
        // Stable sort; the C++ uses std::sort (unstable). Ties among equal
        // x with different deltas are the only place the two can differ.
        intercepts.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap_or(std::cmp::Ordering::Equal));

        let mut winding: i32 = 0;
        let mut run_start: i64 = -1;
        let out_row = r as i32;

        for ix in intercepts.iter() {
            let mut col = ((ix.x - gs.xmin) / dx).floor() as i64;
            if col < 0 {
                col = 0;
            }
            if col >= ncol {
                col = ncol - 1;
            }
            let x_mid = gs.xmin + (col as f64 + 0.5) * dx;
            let left_of_center = ix.x <= x_mid;
            let new_winding = winding + ix.delta;

            if left_of_center {
                if winding != 0 && new_winding == 0 {
                    // Leaving before the cell centre: close the run before this cell.
                    let col_end = col; // exclusive
                    if run_start >= 0 && run_start < col_end {
                        runs.push(GridRun {
                            row: out_row,
                            col_start: run_start as i32,
                            col_end: col_end as i32,
                            id: poly_id,
                        });
                    }
                    run_start = -1;
                } else if winding == 0 && new_winding != 0 {
                    // Entering: this cell is inside.
                    run_start = col;
                }
            } else if winding != 0 && new_winding == 0 {
                // Leaving after the cell centre: include this cell.
                if run_start < 0 {
                    run_start = col;
                }
                let col_end_1 = (col + 1).min(ncol);
                if run_start < col_end_1 {
                    runs.push(GridRun {
                        row: out_row,
                        col_start: run_start as i32,
                        col_end: col_end_1 as i32,
                        id: poly_id,
                    });
                }
                run_start = -1;
            } else if winding == 0 && new_winding != 0 {
                // Entering after the cell centre: start after this cell.
                run_start = col + 1;
            }
            winding = new_winding;
        }

        // Polygon extends beyond the grid's right edge.
        if winding != 0 && run_start >= 0 && run_start < ncol {
            runs.push(GridRun { row: out_row, col_start: run_start as i32, col_end: ncol as i32, id: poly_id });
        }
    }
}
