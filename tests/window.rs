//! `CoverageWindow` and `RowPrefix` against `burn` in coverage mode.

use controlledburn::*;
use proptest::prelude::*;

const N: usize = 24;

/// Dense coverage from `burn` on an N x N unit grid in index space
/// (world y = -row).
fn burn_dense(rings: &[(Vec<f64>, Vec<f64>)]) -> Vec<f64> {
    let grid = GridSpec::new(0.0, -(N as f64), N as f64, 0.0, N as u32, N as u32);
    let rings: Vec<Vec<Coord>> = rings
        .iter()
        .map(|(xs, ys)| {
            let mut r: Vec<Coord> = xs.iter().zip(ys).map(|(&x, &y)| Coord::new(x, -y)).collect();
            r.push(r[0]);
            r
        })
        .collect();
    let res = burn(&[Geometry::Polygon(Polygon::new(rings))], &grid, BurnOptions::coverage()).unwrap();
    let mut out = vec![0.0; N * N];
    for r in &res.runs {
        for c in r.col_start..r.col_end {
            out[r.row as usize * N + c as usize] += 1.0;
        }
    }
    for e in &res.edges {
        out[e.row as usize * N + e.col as usize] += e.fraction as f64;
    }
    out
}

fn window_dense(rings: &[(Vec<f64>, Vec<f64>)]) -> Vec<f64> {
    let mut w = CoverageWindow::default();
    w.reset(0, 0, N, N);
    for (k, (xs, ys)) in rings.iter().enumerate() {
        w.add_ring(xs, ys, if k == 0 { 1.0 } else { -1.0 });
    }
    w.resolve().to_vec()
}

/// A star-shaped (hence simple) ring around (cx, cy).
fn star(cx: f64, cy: f64, radii: &[f64], phase: f64) -> (Vec<f64>, Vec<f64>) {
    let n = radii.len();
    let mut xs = Vec::with_capacity(n);
    let mut ys = Vec::with_capacity(n);
    for (k, r) in radii.iter().enumerate() {
        let a = phase + std::f64::consts::TAU * k as f64 / n as f64;
        xs.push(cx + r * a.cos());
        ys.push(cy + r * a.sin());
    }
    (xs, ys)
}

fn assert_close(a: &[f64], b: &[f64], tol: f64) {
    for (k, (x, y)) in a.iter().zip(b).enumerate() {
        assert!((x - y).abs() <= tol, "cell {k}: {x} vs {y}");
    }
}

#[test]
fn square_on_the_lattice_and_off_it() {
    for (xs, ys) in [
        (vec![2.0, 6.0, 6.0, 2.0], vec![3.0, 3.0, 9.0, 9.0]),
        (vec![2.3, 6.0, 6.0, 2.3], vec![3.0, 3.0, 9.7, 9.7]),
        (vec![2.3, 2.3, 6.1, 6.1], vec![3.2, 9.7, 9.7, 3.2]),
    ] {
        let r = vec![(xs, ys)];
        assert_close(&window_dense(&r), &burn_dense(&r), 1e-6);
    }
}

#[test]
fn a_ring_partly_outside_the_window() {
    // Clipped by the window on three sides; burn clips by its grid.
    let r = vec![(vec![-3.5, 30.2, 12.0], vec![5.5, 2.1, 40.0])];
    assert_close(&window_dense(&r), &burn_dense(&r), 1e-6);
}

#[test]
fn hole() {
    let outer = star(12.0, 12.0, &[9.0, 10.0, 8.5, 9.5, 10.5, 9.0, 8.0], 0.3);
    let hole = star(12.2, 11.7, &[3.0, 4.0, 3.5, 2.5, 3.9], 1.1);
    let r = vec![outer, hole];
    assert_close(&window_dense(&r), &burn_dense(&r), 1e-6);
}

#[test]
fn prefix_integral_equals_coverage_weighted_sum() {
    let field = |i: usize, j: usize| ((i * 7 + j * 3) % 11) as f64 - 4.0;
    let pre = RowPrefix::<2>::from_fn(N, N, |i, j| [field(i, j), 1.0]);
    let (xs, ys) = star(11.0, 12.5, &[6.0, 9.0, 7.5, 10.0, 5.0, 8.0, 9.5, 6.5], 0.2);
    let cov = window_dense(&[(xs.clone(), ys.clone())]);
    let want: f64 = cov.iter().enumerate().map(|(k, c)| c * field(k % N, k / N)).sum();
    let area: f64 = cov.iter().sum();
    let got = pre.integrate_ring(&xs, &ys);
    assert!((got[0] - want).abs() < 1e-9, "{} vs {}", got[0], want);
    assert!((got[1] - area).abs() < 1e-9, "{} vs {}", got[1], area);
    // Reversed ring, same answer.
    let (rx, ry): (Vec<f64>, Vec<f64>) = (xs.iter().rev().cloned().collect(), ys.iter().rev().cloned().collect());
    let back = pre.integrate_ring(&rx, &ry);
    assert!((back[0] - want).abs() < 1e-9);
}

#[test]
fn shared_edges_add_up_to_the_whole() {
    // Two quads sharing the edge b -> c: integrating each ring, or the
    // outer boundary once, gives the same total.
    let pre = RowPrefix::<1>::from_fn(N, N, |i, j| [(i + 2 * j) as f64]);
    let (a, b, c, d, e, f) = ((2.2, 3.1), (9.7, 2.4), (10.4, 11.3), (1.6, 12.2), (18.9, 3.3), (19.5, 12.8));
    let ring = |p: &[(f64, f64)]| -> [f64; 1] {
        let (xs, ys): (Vec<f64>, Vec<f64>) = p.iter().cloned().unzip();
        pre.integrate_ring(&xs, &ys)
    };
    let left = ring(&[a, b, c, d]);
    let right = ring(&[b, e, f, c]);
    let whole = ring(&[a, b, e, f, c, d]);
    assert!((left[0] + right[0] - whole[0]).abs() < 1e-9);
    // Open polylines: the shared edge b -> c, run once each way, cancels,
    // so it can be integrated once and booked +/- to the two rings.
    let line = |p: &[(f64, f64)]| -> f64 {
        let (xs, ys): (Vec<f64>, Vec<f64>) = p.iter().cloned().unzip();
        pre.integrate_polyline(&xs, &ys)[0]
    };
    let bc = line(&[b, c]);
    let left_raw = line(&[c, d, a, b]) + bc;
    let right_raw = line(&[b, e, f, c]) - bc;
    let whole_raw = line(&[a, b, e, f, c, d, a]);
    assert!((left_raw + right_raw - whole_raw).abs() < 1e-9);
    assert!((left_raw.abs() - left[0].abs()).abs() < 1e-9);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(400))]

    #[test]
    fn star_rings_match_burn(
        cx in 4.0f64..20.0, cy in 4.0f64..20.0,
        radii in prop::collection::vec(0.5f64..9.0, 3..24),
        phase in 0.0f64..6.3,
        reverse in any::<bool>(),
    ) {
        let (mut xs, mut ys) = star(cx, cy, &radii, phase);
        if reverse { xs.reverse(); ys.reverse(); }
        let r = vec![(xs, ys)];
        let a = window_dense(&r);
        let b = burn_dense(&r);
        for (x, y) in a.iter().zip(&b) {
            prop_assert!((x - y).abs() <= 2e-6, "{} vs {}", x, y);
        }
    }

    #[test]
    fn lattice_vertices_match_burn(
        pts in prop::collection::vec((0i32..24, 0i32..24), 3..10),
    ) {
        // Integer vertices: every vertex on a cell corner, edges along walls.
        let mut v: Vec<(f64, f64)> = pts.iter().map(|&(x, y)| (x as f64, y as f64)).collect();
        v.dedup();
        // Angular order about an interior point gives a simple ring.
        let n = v.len() as f64;
        let (cx, cy) = (v.iter().map(|p| p.0).sum::<f64>() / n + 1e-3, v.iter().map(|p| p.1).sum::<f64>() / n + 2e-3);
        v.sort_by(|a, b| (a.1 - cy).atan2(a.0 - cx).partial_cmp(&(b.1 - cy).atan2(b.0 - cx)).unwrap());
        v.dedup();
        prop_assume!(v.len() >= 3);
        let (xs, ys): (Vec<f64>, Vec<f64>) = v.into_iter().unzip();
        prop_assume!(ring_area(&xs, &ys).abs() > 1e-9);
        let r = vec![(xs, ys)];
        let a = window_dense(&r);
        let b = burn_dense(&r);
        for (x, y) in a.iter().zip(&b) {
            prop_assert!((x - y).abs() <= 2e-6, "{} vs {}", x, y);
        }
    }
}
