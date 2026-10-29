//! `Measure::CosLatitude`: the ring measure against closed forms and
//! quadrature, and `CoverageWindow` coverage against an independent path
//! (each cell's clip of the ring, measured by `Measure::ring_area`).

use controlledburn::*;
use proptest::prelude::*;

/// Rows of a 0.5 degree grid from 60 S, going south (y down).
fn sphere() -> Measure {
    Measure::cos_latitude((-60.0f64).to_radians(), (-0.5f64).to_radians())
}

fn phi(m: Measure, y: f64) -> f64 {
    match m {
        Measure::CosLatitude { phi0, dphi } => phi0 + dphi * y,
        _ => unreachable!(),
    }
}

/// A star-shaped (hence simple) ring around (cx, cy).
fn star(cx: f64, cy: f64, radii: &[f64], phase: f64) -> (Vec<f64>, Vec<f64>) {
    let n = radii.len();
    let (mut xs, mut ys) = (Vec::with_capacity(n), Vec::with_capacity(n));
    for (k, r) in radii.iter().enumerate() {
        let a = phase + std::f64::consts::TAU * k as f64 / n as f64;
        xs.push(cx + r * a.cos());
        ys.push(cy + r * a.sin());
    }
    (xs, ys)
}

/// Sutherland-Hodgman clip of a ring to the cell `[i, i + 1] x [j, j + 1]`.
fn clip(xs: &[f64], ys: &[f64], i: f64, j: f64) -> (Vec<f64>, Vec<f64>) {
    let mut pts: Vec<(f64, f64)> = xs.iter().copied().zip(ys.iter().copied()).collect();
    // inside(p) >= 0 for each half-plane; cross at t along p -> q
    let planes: [(f64, f64, f64); 4] = [(1.0, 0.0, -i), (-1.0, 0.0, i + 1.0), (0.0, 1.0, -j), (0.0, -1.0, j + 1.0)];
    for (a, b, c) in planes {
        let f = |p: (f64, f64)| a * p.0 + b * p.1 + c;
        let mut out = Vec::new();
        for k in 0..pts.len() {
            let (p, q) = (pts[k], pts[(k + 1) % pts.len()]);
            let (fp, fq) = (f(p), f(q));
            if fp >= 0.0 {
                out.push(p);
            }
            if (fp >= 0.0) != (fq >= 0.0) {
                let t = fp / (fp - fq);
                out.push((p.0 + t * (q.0 - p.0), p.1 + t * (q.1 - p.1)));
            }
        }
        pts = out;
        if pts.is_empty() {
            break;
        }
    }
    pts.into_iter().unzip()
}

#[test]
fn a_lonlat_rectangle_has_the_spherical_cell_area() {
    let m = sphere();
    let (x0, x1, y0, y1) = (2.0, 9.5, 3.25, 17.0);
    let (xs, ys) = ([x0, x1, x1, x0], [y0, y0, y1, y1]);
    let Measure::CosLatitude { dphi, .. } = m else { unreachable!() };
    let want = (x1 - x0) * (phi(m, y1).sin() - phi(m, y0).sin()) / dphi;
    let got = m.ring_area(&xs, &ys);
    assert!((got - want).abs() < 1e-13 * want.abs(), "{got} vs {want}");
    // same sign as the planar ring_area
    assert_eq!(got.signum(), ring_area(&xs, &ys).signum());
    // and a whole cell is cell_area
    let (xs, ys) = ([4.0, 5.0, 5.0, 4.0], [7.0, 7.0, 8.0, 8.0]);
    assert!((m.ring_area(&xs, &ys) - m.cell_area(7)).abs() < 1e-15);
}

#[test]
fn ring_measure_matches_quadrature() {
    // A triangle with slanted edges, integrated on a fine midpoint grid.
    let m = sphere();
    let (xs, ys) = ([1.0, 30.0, 9.0], [2.0, 11.0, 60.0]);
    let inside = |x: f64, y: f64| {
        let s = |k: usize| {
            let k1 = (k + 1) % 3;
            (xs[k1] - xs[k]) * (y - ys[k]) - (ys[k1] - ys[k]) * (x - xs[k])
        };
        let (a, b, c) = (s(0), s(1), s(2));
        (a >= 0.0 && b >= 0.0 && c >= 0.0) || (a <= 0.0 && b <= 0.0 && c <= 0.0)
    };
    let n = 3000;
    let (h_x, h_y) = (30.0 / n as f64, 60.0 / n as f64);
    let mut q = 0.0;
    for iy in 0..n {
        let y = (iy as f64 + 0.5) * h_y;
        let w = phi(m, y).cos();
        for ix in 0..n {
            let x = (ix as f64 + 0.5) * h_x;
            if inside(x, y) {
                q += w;
            }
        }
    }
    q *= h_x * h_y;
    let got = m.ring_area(&xs, &ys).abs();
    assert!((got - q).abs() < 2e-3 * q, "{got} vs quadrature {q}");
}

#[test]
fn small_dphi_tends_to_planar_times_cos() {
    let p0 = 0.7f64;
    let m = Measure::cos_latitude(p0, 1e-9);
    let (xs, ys) = star(5.0, 5.0, &[3.0, 2.0, 3.5, 1.5, 2.5], 0.3);
    let r = m.ring_area(&xs, &ys) / (ring_area(&xs, &ys) * p0.cos());
    assert!((r - 1.0).abs() < 1e-7, "{r}");
}

/// Coverage per cell, from the window, against the clipped ring's measure
/// over the cell's.
fn check_against_clip(m: Measure, xs: &[f64], ys: &[f64], tol: f64) {
    let mut w = CoverageWindow::with_measure(m);
    assert!(w.reset_to_ring(xs, ys));
    w.add_ring(xs, ys, 1.0);
    let cov = w.resolve().to_vec();
    let (c0, r0, nc) = (w.col0(), w.row0(), w.ncol());
    let mut total = 0.0;
    for (k, &c) in cov.iter().enumerate() {
        let (i, j) = ((c0 + (k % nc) as i64) as f64, (r0 + (k / nc) as i64) as f64);
        let (cx, cy) = clip(xs, ys, i, j);
        let want = if cx.len() < 3 { 0.0 } else { m.ring_area(&cx, &cy).abs() / m.cell_area(j as i64) };
        assert!((c - want).abs() <= tol, "cell ({i}, {j}): {c} vs {want}");
        total += c * m.cell_area(j as i64);
    }
    let whole = m.ring_area(xs, ys).abs();
    assert!((total - whole).abs() <= tol * whole, "total {total} vs {whole}");
}

#[test]
fn coverage_matches_clipped_measure_on_a_square_and_a_star() {
    let m = sphere();
    check_against_clip(m, &[0.5, 2.5, 2.5, 0.5], &[0.5, 0.5, 2.5, 2.5], 1e-13);
    let (xs, ys) = star(12.3, 40.7, &[6.0, 3.0, 7.5, 2.0, 5.0, 6.5, 4.0], 0.2);
    check_against_clip(m, &xs, &ys, 1e-12);
    // the other way round, and a hole cancels its part
    let (rx, ry): (Vec<f64>, Vec<f64>) = (xs.iter().rev().copied().collect(), ys.iter().rev().copied().collect());
    check_against_clip(m, &rx, &ry, 1e-12);
}

#[test]
fn planar_measure_is_unchanged() {
    let (xs, ys) = star(7.3, 6.1, &[4.0, 2.0, 5.0, 1.5, 3.0], 0.4);
    let mut a = CoverageWindow::default();
    let mut b = CoverageWindow::with_measure(Measure::Planar);
    for w in [&mut a, &mut b] {
        w.reset(0, 0, 16, 16);
        w.add_ring(&xs, &ys, 1.0);
        w.resolve();
    }
    assert_eq!(a.coverage(), b.coverage());
    assert_eq!(Measure::Planar.ring_area(&xs, &ys), ring_area(&xs, &ys));
}

proptest! {
    #[test]
    fn star_rings_match_clipped_measure(
        cx in 4.0f64..40.0, cy in 4.0f64..100.0,
        radii in prop::collection::vec(0.3f64..3.5, 3..12),
        phase in 0.0f64..6.3,
        phi0 in -80.0f64..80.0, dphi in -2.0f64..2.0,
    ) {
        prop_assume!(dphi.abs() > 1e-3);
        // keep latitudes within the poles over the ring
        let m = Measure::cos_latitude(phi0.to_radians(), dphi.to_radians() / 10.0);
        let (lo, hi) = (phi(m, cy - 4.0).to_degrees(), phi(m, cy + 4.0).to_degrees());
        prop_assume!(lo.abs() < 89.0 && hi.abs() < 89.0);
        let (xs, ys) = star(cx, cy, &radii, phase);
        check_against_clip(m, &xs, &ys, 1e-10);
    }
}

#[test]
fn long_edges_equal_their_subdivision() {
    // Edges spanning 160 degrees of latitude take the direct branch of the
    // edge integral; split into short pieces they take the series.
    let m = Measure::cos_latitude(80f64.to_radians(), (-1f64).to_radians());
    let (xs, ys) = ([0.0, 30.0, 7.0], [0.0, 160.0, 90.0]);
    let (mut fx, mut fy) = (Vec::new(), Vec::new());
    let n = 4000;
    for k in 0..3 {
        let k1 = (k + 1) % 3;
        for s in 0..n {
            let t = s as f64 / n as f64;
            fx.push(xs[k] + t * (xs[k1] - xs[k]));
            fy.push(ys[k] + t * (ys[k1] - ys[k]));
        }
    }
    let (a, b) = (m.ring_area(&xs, &ys), m.ring_area(&fx, &fy));
    assert!((a - b).abs() < 1e-11 * a.abs(), "{a} vs {b}");
}
