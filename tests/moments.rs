//! `CoverageWindow`'s first moments against an independent path: each
//! cell's clip of the ring, its moments by the polygon formula (planar) or
//! by Green's theorem with Gauss-Legendre quadrature along each edge
//! (`Measure::CosLatitude`).

use controlledburn::*;
use proptest::prelude::*;

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
fn clip(xs: &[f64], ys: &[f64], i: f64, j: f64) -> Vec<(f64, f64)> {
    let mut pts: Vec<(f64, f64)> = xs.iter().copied().zip(ys.iter().copied()).collect();
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
    pts
}

/// 16-point Gauss-Legendre nodes and weights on [-1, 1] (positive half).
const GL: [(f64, f64); 8] = [
    (0.0950125098376374, 0.1894506104550685),
    (0.2816035507792589, 0.1826034150449236),
    (0.4580167776572274, 0.1691565193950025),
    (0.6178762444026438, 0.1495959888165767),
    (0.755404408355003, 0.1246289712555339),
    (0.8656312023878318, 0.0951585116824928),
    (0.9445750230732326, 0.0622535239386479),
    (0.9894009349916499, 0.0271524594117541),
];

/// Measure, and first moments about `(i, j)`, of a polygon under `w(y)`,
/// by Green's theorem with quadrature along each edge: `int w dy x`,
/// `int (x - i)^2 / 2 w dy`, `int (x - i) (y - j) w dy`. Signed by the
/// polygon's orientation.
fn moments_by_quadrature(pts: &[(f64, f64)], i: f64, j: f64, w: &dyn Fn(f64) -> f64) -> [f64; 3] {
    let mut out = [0.0; 3];
    let n = pts.len();
    for k in 0..n {
        let (p, q) = (pts[k], pts[(k + 1) % n]);
        let (dx, dy) = (q.0 - p.0, q.1 - p.1);
        for &(t, wt) in &GL {
            for t in [t, -t] {
                let s = 0.5 * (1.0 + t);
                let (x, y) = (p.0 + s * dx, p.1 + s * dy);
                let f = 0.5 * wt * dy * w(y);
                out[0] += (x - i) * f;
                out[1] += 0.5 * (x - i) * (x - i) * f;
                out[2] += (x - i) * (y - j) * f;
            }
        }
    }
    out
}

/// Planar measure and moments about `(i, j)` by the polygon formulas.
fn moments_planar(pts: &[(f64, f64)], i: f64, j: f64) -> [f64; 3] {
    let n = pts.len();
    let mut out = [0.0; 3];
    for k in 0..n {
        let (p, q) = (pts[k], pts[(k + 1) % n]);
        let (x0, y0, x1, y1) = (p.0 - i, p.1 - j, q.0 - i, q.1 - j);
        // with y down, positive area is clockwise as drawn: x dy
        let c = x0 * y1 - x1 * y0;
        out[0] += 0.5 * c;
        out[1] += (x0 + x1) * c / 6.0;
        out[2] += (y0 + y1) * c / 6.0;
    }
    out
}

/// Each cell's coverage and moments from the window, against the clip's.
fn check(m: Measure, rings: &[(Vec<f64>, Vec<f64>)], tol: f64) {
    let mut w = CoverageWindow::with_measure(m);
    w.set_moments(true);
    let (xs, ys) = &rings[0];
    assert!(w.reset_to_ring(xs, ys));
    for (k, (xs, ys)) in rings.iter().enumerate() {
        w.add_ring(xs, ys, if k == 0 { 1.0 } else { -1.0 });
    }
    w.resolve();
    let (c0, r0, nc) = (w.col0(), w.row0(), w.ncol());
    let weight = |y: f64| match m {
        Measure::CosLatitude { phi0, dphi } => (phi0 + dphi * y).cos(),
        _ => 1.0,
    };
    for k in 0..w.coverage().len() {
        let (i, j) = ((c0 + (k % nc) as i64) as f64, (r0 + (k / nc) as i64) as f64);
        let mut want = [0.0; 3];
        for (r, (xs, ys)) in rings.iter().enumerate() {
            let pts = clip(xs, ys, i, j);
            if pts.len() < 3 {
                continue;
            }
            let v = match m {
                Measure::Planar => moments_planar(&pts, i, j),
                _ => moments_by_quadrature(&pts, i, j, &weight),
            };
            // exterior positive, holes negative, whichever way each runs
            let s = if r == 0 { 1.0 } else { -1.0 } * v[0].signum();
            for t in 0..3 {
                want[t] += s * v[t];
            }
        }
        let cell = m.cell_area(j as i64);
        let got = [w.coverage()[k], w.moment_x()[k], w.moment_y()[k]];
        for t in 0..3 {
            let e = (got[t] - want[t] / cell).abs();
            assert!(e <= tol, "cell ({i}, {j}) term {t}: {} vs {} ({e:e})", got[t], want[t] / cell);
        }
    }
}

fn sphere() -> Measure {
    Measure::cos_latitude((-60.0f64).to_radians(), (-0.5f64).to_radians())
}

#[test]
fn a_whole_cell_has_its_centre_as_centroid() {
    for m in [Measure::Planar, sphere()] {
        let mut w = CoverageWindow::with_measure(m);
        w.set_moments(true);
        w.reset(0, 0, 5, 5);
        w.add_ring(&[1.0, 4.0, 4.0, 1.0], &[1.0, 1.0, 4.0, 4.0], 1.0);
        w.resolve();
        let k = 2 * 5 + 2;
        assert!((w.coverage()[k] - 1.0).abs() < 1e-15);
        assert!((w.moment_x()[k] - 0.5).abs() < 1e-15, "{}", w.moment_x()[k]);
        // the planar centre is the middle; on the sphere it moves to the
        // equator's side (rows go south from 60 S, so up, to smaller y)
        let my = w.moment_y()[k];
        match m {
            Measure::Planar => assert!((my - 0.5).abs() < 1e-15),
            _ => assert!(my < 0.5 && my > 0.49, "{my}"),
        }
    }
}

#[test]
fn moments_match_the_clip_on_a_square_a_star_and_a_hole() {
    for (m, tol) in [(Measure::Planar, 1e-13), (sphere(), 1e-12)] {
        check(m, &[(vec![0.5, 2.5, 2.5, 0.5], vec![0.5, 0.5, 2.5, 2.5])], tol);
        let star_ = star(12.3, 40.7, &[6.0, 3.0, 7.5, 2.0, 5.0, 6.5, 4.0], 0.2);
        check(m, std::slice::from_ref(&star_), tol);
        let rev: (Vec<f64>, Vec<f64>) =
            (star_.0.iter().rev().copied().collect(), star_.1.iter().rev().copied().collect());
        check(m, &[rev], tol);
        let hole = star(12.0, 40.0, &[1.5, 1.0, 1.7, 0.8, 1.2], 0.9);
        check(m, &[star_, hole], tol);
    }
}

#[test]
fn moments_leave_coverage_as_it_was() {
    let (xs, ys) = star(7.3, 6.1, &[4.0, 2.0, 5.0, 1.5, 3.0], 0.4);
    for m in [Measure::Planar, sphere()] {
        for weight in [1.0, -1.0, 0.37] {
            let mut a = CoverageWindow::with_measure(m);
            let mut b = CoverageWindow::with_measure(m);
            b.set_moments(true);
            for w in [&mut a, &mut b] {
                w.reset(0, 0, 16, 16);
                w.add_ring(&xs, &ys, weight);
                w.resolve();
            }
            assert_eq!(a.coverage(), b.coverage());
            assert!(a.moment_x().is_empty() && a.moment_y().is_empty());
            assert_eq!(b.moment_x().len(), 256);
        }
    }
}

#[test]
fn a_ring_partly_outside_the_window() {
    // pieces left and right of the window clamp to its walls; the cells
    // inside keep their moments
    let (xs, ys) = star(5.2, 5.4, &[4.0, 2.5, 4.5, 3.0, 3.5, 2.0], 0.1);
    for m in [Measure::Planar, sphere()] {
        let mut full = CoverageWindow::with_measure(m);
        let mut part = CoverageWindow::with_measure(m);
        for w in [&mut full, &mut part] {
            w.set_moments(true);
        }
        full.reset(0, 0, 12, 12);
        part.reset(3, 2, 4, 6);
        for w in [&mut full, &mut part] {
            w.add_ring(&xs, &ys, 1.0);
            w.resolve();
        }
        for r in 0..6 {
            for c in 0..4 {
                let (k, kf) = (r * 4 + c, (r + 2) * 12 + c + 3);
                for (a, b) in [(part.moment_x()[k], full.moment_x()[kf]), (part.moment_y()[k], full.moment_y()[kf])] {
                    assert!((a - b).abs() < 1e-13, "cell ({}, {}): {a} vs {b}", c + 3, r + 2);
                }
            }
        }
    }
}

#[test]
fn moments_add_up_to_the_rings() {
    // over the cells, i * coverage + moment_x (times the cell's measure) is
    // the ring's own first moment
    let (xs, ys) = star(9.6, 31.2, &[6.0, 3.0, 7.0, 2.5, 5.5], 0.7);
    for m in [Measure::Planar, sphere()] {
        let mut w = CoverageWindow::with_measure(m);
        w.set_moments(true);
        assert!(w.reset_to_ring(&xs, &ys));
        w.add_ring(&xs, &ys, 1.0);
        w.resolve();
        let (mut mx, mut my) = (0.0, 0.0);
        for k in 0..w.coverage().len() {
            let (i, j) = ((w.col0() + (k % w.ncol()) as i64) as f64, (w.row0() + (k / w.ncol()) as i64) as f64);
            let cell = m.cell_area(j as i64);
            mx += (i * w.coverage()[k] + w.moment_x()[k]) * cell;
            my += (j * w.coverage()[k] + w.moment_y()[k]) * cell;
        }
        let pts: Vec<(f64, f64)> = xs.iter().copied().zip(ys.iter().copied()).collect();
        let weight = |y: f64| match m {
            Measure::CosLatitude { phi0, dphi } => (phi0 + dphi * y).cos(),
            _ => 1.0,
        };
        let v = moments_by_quadrature(&pts, 0.0, 0.0, &weight);
        let s = v[0].signum();
        assert!((mx - s * v[1]).abs() < 1e-10 * v[1].abs(), "{mx} vs {}", s * v[1]);
        assert!((my - s * v[2]).abs() < 1e-10 * v[2].abs(), "{my} vs {}", s * v[2]);
    }
}

#[test]
fn tall_rows_take_the_closed_forms() {
    // rows of 1.2 radians: a piece spanning most of a row takes the direct
    // branch of the series, not the expansion
    let m = Measure::cos_latitude(-1.3, 1.2);
    check(m, &[(vec![0.3, 2.7, 1.9, 0.1], vec![0.05, 0.2, 2.05, 1.9])], 1e-11);
}

proptest! {
    #[test]
    fn star_rings_match_the_clip(
        cx in 4.0f64..40.0, cy in 4.0f64..100.0,
        radii in prop::collection::vec(0.3f64..3.5, 3..12),
        phase in 0.0f64..6.3,
        phi0 in -80.0f64..80.0, dphi in -2.0f64..2.0,
        planar in any::<bool>(),
    ) {
        let m = if planar {
            Measure::Planar
        } else {
            prop_assume!(dphi.abs() > 1e-3);
            let m = Measure::cos_latitude(phi0.to_radians(), dphi.to_radians() / 10.0);
            let Measure::CosLatitude { phi0, dphi } = m else { unreachable!() };
            let (lo, hi) = ((phi0 + dphi * (cy - 4.0)).to_degrees(), (phi0 + dphi * (cy + 4.0)).to_degrees());
            prop_assume!(lo.abs() < 89.0 && hi.abs() < 89.0);
            m
        };
        check(m, &[star(cx, cy, &radii, phase)], 1e-10);
    }
}
