//! Approx-mode polygon rasterization: the cell-centre rule (fasterize
//! semantics) via a direct edge/row-centre intersection sweep. No walker,
//! no coverage math. O(edges * rows touched).
//!
//! A cell is inside iff the winding number at its centre is non-zero,
//! with a left-inclusive convention: an edge crossing exactly at the cell
//! centre x counts as left of centre. Which side of a centre an edge
//! passes is decided exactly (an orientation test on the edge's
//! endpoints), so the rule holds for the centre and edges as given.

use crate::ee::bbox::BBox;
use crate::geometry::{is_ccw, Polygon};
use crate::grid::GridSpec;
use crate::output::GridRun;
use crate::polygon::ring_envelope;

#[derive(Clone, Copy)]
struct Intercept {
    x: f64,
    /// `4 * col + 2 * side + up`: the column (clamped to the grid), `side`
    /// 0 if the crossing is left of or on that column's centre (decided
    /// exactly) and 1 if right of it, and `up` 1 for a winding delta of +1,
    /// 0 for -1. Sorting by `key >> 1` orders by column, then side.
    key: i64,
}

/// `a + b` as the rounded sum and its exact error.
#[inline]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    (s, (a - (s - bb)) + (b - bb))
}

/// Sign of `(bx - ax)(cy - ay) - (cx - ax)(by - ay)`, exactly: positive
/// when `c` is left of the line `a -> b`. The determinant is expanded into
/// six products of input coordinates, each split exactly with an FMA, and
/// the twelve terms are summed as a non-overlapping expansion (Shewchuk's
/// grow-expansion), whose sign is that of its largest component. A float
/// filter answers the common case.
pub(crate) fn orient2d_sign(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> i32 {
    let l = (b.0 - a.0) * (c.1 - a.1);
    let r = (c.0 - a.0) * (b.1 - a.1);
    let det = l - r;
    // Shewchuk's ccwerrboundA, padded.
    let bound = 4.0e-16 * (l.abs() + r.abs());
    if det > bound {
        return 1;
    }
    if det < -bound {
        return -1;
    }
    let prods =
        [(b.0, c.1, 1.0), (b.0, a.1, -1.0), (a.0, c.1, -1.0), (c.0, b.1, -1.0), (c.0, a.1, 1.0), (a.0, b.1, 1.0)];
    let mut e: [f64; 12] = [0.0; 12];
    let mut n = 0;
    for &(p, q, sgn) in &prods {
        let hi = p * q;
        let lo = p.mul_add(q, -hi);
        for t in [sgn * lo, sgn * hi] {
            // grow-expansion: add t to e[..n], keeping components in
            // increasing magnitude
            let mut q = t;
            for ek in e.iter_mut().take(n) {
                let (s, err) = two_sum(q, *ek);
                *ek = err;
                q = s;
            }
            e[n] = q;
            n += 1;
        }
    }
    match e[..n].iter().rev().find(|&&v| v != 0.0) {
        Some(&v) if v > 0.0 => 1,
        Some(_) => -1,
        None => 0,
    }
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
                // Crossing direction from the edge's own endpoints, not from
                // y0 against the recomputed y_mid: the row-range test above
                // admits an edge whose top endpoint is on the centre line
                // (within 1e-10), but that endpoint can sit an ulp below
                // y_mid as recomputed from ymax and dy, and `y0 >= y_mid`
                // then tags a downward edge as upward. Its intercept gets the
                // wrong sign, the winding never returns to zero, and the row
                // fills to the grid edge.
                let delta = if y1 < y0 { -1 } else { 1 } * wf;
                // Left of or on the centre (the left-inclusive rule),
                // decided exactly from the edge's endpoints taken bottom to
                // top, not from the rounded intercept: two edges through
                // nearly the same point (a sliver, or a ring doubling back
                // along itself) give intercepts that round either side of
                // a centre they both miss, and the cell burned. The rounded
                // intercept is within about 8 eps (|x0| + |x1|) of the exact
                // one; outside a generous multiple of that, the float
                // comparison already gives the exact answer.
                let col = (((x_int - gs.xmin) / dx).floor() as i64).clamp(0, ncol - 1);
                let x_mid = gs.xmin + (col as f64 + 0.5) * dx;
                let left = if (x_int - x_mid).abs() > 1e-13 * (x0.abs() + x1.abs()) {
                    x_int < x_mid
                } else {
                    let (lo, hi) = if y0 < y1 { ((x0, y0), (x1, y1)) } else { ((x1, y1), (x0, y0)) };
                    orient2d_sign(lo, hi, (x_mid, y_mid)) <= 0
                };
                let key = 4 * col + 2 * i64::from(!left) + i64::from(delta > 0);
                row_intercepts[(r - row_lo) as usize].push(Intercept { x: x_int, key });
            }
        }
    }

    for r in row_lo..=row_hi {
        let intercepts = &mut row_intercepts[(r - row_lo) as usize];
        if intercepts.is_empty() {
            continue;
        }
        // Ordered by column, then side of its centre, then x: only the side
        // of each centre matters, and this keeps the order consistent with
        // the exact side test where rounded intercepts tie or cross. Sorted
        // by x, then a stable insertion pass on (column, side), which the x
        // order already follows except within rounding of a centre: the
        // same order as the C++'s std::stable_sort on (column, side, x), at
        // the cost of one sort by x.
        intercepts.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap_or(std::cmp::Ordering::Equal));
        for i in 1..intercepts.len() {
            let mut j = i;
            while j > 0 && (intercepts[j - 1].key >> 1) > (intercepts[j].key >> 1) {
                intercepts.swap(j - 1, j);
                j -= 1;
            }
        }

        let mut winding: i32 = 0;
        let mut run_start: i64 = -1;
        let out_row = r as i32;

        for ix in intercepts.iter() {
            let col = ix.key >> 2;
            let left_of_center = ix.key & 2 == 0;
            let delta = if ix.key & 1 == 1 { 1 } else { -1 };
            let new_winding = winding + delta;

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

#[cfg(test)]
mod tests {
    use super::orient2d_sign;

    #[test]
    fn orient2d_is_exact_near_collinear() {
        // The tie from controlledburn-rs#6: the centre (1000.85, -2000.25)
        // against an edge whose rounded intercept is exactly 1000.85 but
        // whose exact intercept is 4.9e-14 to the right of it.
        let lo = (1000.4, -2000.7);
        let hi = (1001.45, -1999.65);
        assert_eq!(orient2d_sign(lo, hi, (1000.85, -2000.25)), 1);
        // exactly collinear, and either side by one ulp
        assert_eq!(orient2d_sign((0.0, 0.0), (1.0, 1.0), (0.5, 0.5)), 0);
        assert_eq!(orient2d_sign((0.0, 0.0), (1.0, 1.0), (0.5, f64::from_bits(0.5f64.to_bits() + 1))), 1);
        assert_eq!(orient2d_sign((0.0, 0.0), (1.0, 1.0), (0.5, f64::from_bits(0.5f64.to_bits() - 1))), -1);
        // Shewchuk's classic near-degenerate grid: points 0.5 + i * 2^-53
        // against the line through (12, 12) and (24, 24); the sign is that
        // of y - x, exactly.
        let ulp = 2f64.powi(-53);
        for i in 0..64 {
            for j in 0..64 {
                let c = (0.5 + i as f64 * ulp, 0.5 + j as f64 * ulp);
                let want = (j - i).signum();
                assert_eq!(orient2d_sign((12.0, 12.0), (24.0, 24.0), c), want, "{i} {j}");
            }
        }
    }
}
