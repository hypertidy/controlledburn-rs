//! Approx-mode polygon rasterization: the cell-centre rule (fasterize
//! semantics) via a direct edge/row-centre intersection sweep. No walker,
//! no coverage math. O(edges * rows touched).
//!
//! A cell is inside iff the winding number at its centre is non-zero,
//! with a left-inclusive convention: an edge crossing exactly at the cell
//! centre x counts as left of centre.

use crate::ee::bbox::BBox;
use crate::geometry::{is_ccw, Polygon};
use crate::grid::Grid;
use crate::output::GridRun;
use crate::polygon::ring_envelope;

#[derive(Clone, Copy)]
struct Intercept {
    x: f64,
    delta: i32,
}

pub(crate) fn process_polygon_approx(poly: &Polygon, grid: &Grid<0>, poly_id: i32, runs: &mut Vec<GridRun>) {
    let Some(exterior) = poly.rings.first() else { return };
    if exterior.len() < 4 {
        return;
    }
    let (xa, ya) = (grid.x_axis(), grid.y_axis());
    let ncol = xa.n() as i64;
    let nrow = ya.n() as i64;

    let env = ring_envelope(exterior);
    let grid_box: BBox = *grid.extent();
    if !env.intersects(&grid_box) {
        return;
    }

    // Row range the polygon can touch (0-based, row 0 at the top).
    let row_lo = (ya.units(env.ymax).floor() as i64).max(0);
    let row_hi = (ya.units(env.ymin).floor() as i64).min(nrow - 1);
    if row_lo > row_hi {
        return;
    }
    let ya_units = |v: f64| ya.units(v);
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
            let (ya_e, yb) = if y0 > y1 { (y1, y0) } else { (y0, y1) };

            // Rows whose y_mid lies in (ya, yb]: top-inclusive,
            // bottom-exclusive, matching the walker's half-open crossing
            // convention. y_mid(r) is the centre of row r, at units r + 0.5.
            let e_row_lo = ((ya_units(yb) - 0.5 - 1e-10).ceil() as i64).max(row_lo);
            let e_row_hi = (((ya_units(ya_e) - 0.5 - 1e-10).ceil() as i64) - 1).min(row_hi);
            if e_row_lo > e_row_hi {
                continue;
            }

            let inv_dy_edge = 1.0 / (y1 - y0);
            for r in e_row_lo..=e_row_hi {
                let y_mid = ya.centre(r);
                let t = (y_mid - y0) * inv_dy_edge;
                let x_int = x0 + t * (x1 - x0);
                // Crossing direction from the edge's own endpoints, not from
                // y0 against the recomputed y_mid: the row-range test above
                // admits an edge whose top endpoint is on the centre line
                // (within 1e-10), but that endpoint can sit an ulp below
                // y_mid as recomputed from ymax and dy, and `y0 >= y_mid`
                // then tags a downward edge as upward. Its intercept gets the
                // wrong sign, the winding never returns to zero, and the row
                // fills to the grid edge.
                let delta = if y1 < y0 { -1 } else { 1 } * wf;
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
            let mut col = xa.units(ix.x).floor() as i64;
            if col < 0 {
                col = 0;
            }
            if col >= ncol {
                col = ncol - 1;
            }
            let x_mid = xa.centre(col);
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
