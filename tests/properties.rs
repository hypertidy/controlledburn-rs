//! Property-based tests: lattice-biased random rings on non-dyadic grids,
//! checked against invariants that must hold for every input.
//!
//! The generator deliberately steers vertices onto the places where parsed
//! input and computed cell edges disagree by an ulp, and where traversal
//! logic has to make a tie-breaking decision: cell corners, cell walls,
//! row and column centre lines, the grid boundary, and just outside it.
//! Uniform random vertices essentially never hit these seams; both bugs
//! fixed on this branch (walker extrapolation off a near-edge vertex,
//! winding loss on a row-centre column-wall crossing) needed one.
//!
//! Properties:
//!   1. coverage conserves area: sum(runs) + sum(edge fractions), times the
//!      cell area, equals the planar area of the ring clipped to the grid.
//!   2. the tables are well formed: fractions in (0, 1], runs inside the
//!      grid with col_start < col_end, no notes.
//!   3. the per-cell cover matrix is invariant under ring orientation and
//!      vertex start rotation.
//!   4. the per-cell cover matrix matches an exact per-cell oracle (the
//!      ring clipped to each cell, shoelace) on small grids; this localises
//!      a failure to a cell where property 1 only reports a total.
//!   5. a cell fully covered in coverage mode is in approx mode's runs and a
//!      cell untouched in coverage mode is not.
//!   6. burning at 2x resolution and averaging 2x2 blocks equals burning at
//!      1x (right in total but wrong in placement is still a failure).
//!
//! Run a long campaign with
//!   PROPTEST_CASES=20000 cargo test --test properties --release

mod common;

use common::*;
use controlledburn::*;
use proptest::prelude::*;

// ---------------------------------------------------------------- grids --

#[derive(Clone, Copy, Debug)]
struct Ext {
    xmin: f64,
    ymin: f64,
    xmax: f64,
    ymax: f64,
}

#[derive(Clone, Debug)]
struct G {
    spec: GridSpec,
}

impl G {
    fn ext(&self) -> Ext {
        let e = self.spec.extent();
        Ext { xmin: e.xmin, ymin: e.ymin, xmax: e.xmax, ymax: e.ymax }
    }
    fn res(&self) -> f64 {
        self.spec.dx()
    }
    fn ncol(&self) -> usize {
        ((self.spec.extent().xmax - self.spec.extent().xmin) / self.spec.dx()).round() as usize
    }
    fn nrow(&self) -> usize {
        ((self.spec.extent().ymax - self.spec.extent().ymin) / self.spec.dy()).round() as usize
    }
}

/// Non-dyadic resolutions on assorted origins, so that `origin + k * res`
/// rarely equals the decimal a vertex was given as.
fn grid_up_to(max_dim: u32) -> impl Strategy<Value = G> {
    let res = prop::sample::select(vec![0.1, 0.05, 0.3, 1.0 / 3.0, 1.0 / 7.0, 0.7, 1.0, 2.5]);
    let origin = prop::sample::select(vec![(0.0, 0.0), (-180.0, -90.0), (0.3, 0.3), (-1.7, 2.9), (1000.1, -2000.7)]);
    (res, origin, 4u32..max_dim, 4u32..max_dim).prop_map(|(res, (x0, y0), ncol, nrow)| G {
        spec: GridSpec::new(x0, y0, x0 + ncol as f64 * res, y0 + nrow as f64 * res, ncol, nrow),
    })
}

fn grid() -> impl Strategy<Value = G> {
    grid_up_to(40)
}

// ---------------------------------------------------------------- rings --

/// Where a vertex is steered to, after a uniform draw inside (or a little
/// outside) the grid.
#[derive(Clone, Copy, Debug)]
enum Snap {
    None,
    Corner,
    XWall,
    YWall,
    RowCentre,
    ColCentre,
    Centre,
    GridEdge,
    Outside,
}

fn snap_kind() -> impl Strategy<Value = Snap> {
    prop_oneof![
        3 => Just(Snap::None),
        3 => Just(Snap::Corner),
        2 => Just(Snap::XWall),
        2 => Just(Snap::YWall),
        2 => Just(Snap::RowCentre),
        2 => Just(Snap::ColCentre),
        1 => Just(Snap::Centre),
        1 => Just(Snap::GridEdge),
        1 => Just(Snap::Outside),
    ]
}

fn apply_snap(g: &G, x: f64, y: f64, s: Snap) -> (f64, f64) {
    let e = g.spec.extent();
    let r = g.res();
    // Compute lattice positions as a decimal a user might type: round the
    // multiplier, then recompute from the origin. That reproduces the
    // "0.3 vs 3 * 0.1" disagreement rather than copying the grid's value.
    let kx = ((x - e.xmin) / r).round();
    let ky = ((e.ymax - y) / r).round();
    let xl = e.xmin + kx * r;
    let yl = e.ymax - ky * r;
    let xc = e.xmin + (kx + 0.5) * r;
    let yc = e.ymax - (ky + 0.5) * r;
    match s {
        Snap::None => (x, y),
        Snap::Corner => (xl, yl),
        Snap::XWall => (xl, y),
        Snap::YWall => (x, yl),
        Snap::RowCentre => (x, yc),
        Snap::ColCentre => (xc, y),
        Snap::Centre => (xc, yc),
        Snap::GridEdge => {
            // snap to the nearest grid boundary in one axis
            let dl = (x - e.xmin).abs();
            let dr = (x - e.xmax).abs();
            let db = (y - e.ymin).abs();
            let dt = (y - e.ymax).abs();
            let m = dl.min(dr).min(db).min(dt);
            if m == dl {
                (e.xmin, y)
            } else if m == dr {
                (e.xmax, y)
            } else if m == db {
                (x, e.ymin)
            } else {
                (x, e.ymax)
            }
        }
        Snap::Outside => {
            let w = e.xmax - e.xmin;
            let h = e.ymax - e.ymin;
            // push past the nearest boundary by up to one grid width
            if (x - e.xmin) < w / 2.0 {
                (e.xmin - (x - e.xmin).abs() - 0.37 * r, y)
            } else {
                (x, e.ymax + (y - e.ymin).abs() * 0.5 + 0.37 * h / g.nrow() as f64)
            }
        }
    }
}

/// A simple (non-self-intersecting) ring: vertices sorted by angle around
/// a centre, radii random, then each vertex snapped. Snapping can in
/// principle create a self-intersection for very spiky rings; those are
/// filtered with a segment test so every ring is a valid polygon.
fn ring_on(g: G) -> impl Strategy<Value = (G, Vec<(f64, f64)>)> {
    let e = g.spec.extent();
    let w = e.xmax - e.xmin;
    let h = e.ymax - e.ymin;
    let n = 3usize..12;
    let verts = n.prop_flat_map(move |n| {
        ((0.05f64..0.95, 0.05f64..0.95), prop::collection::vec((0.0f64..1.0, 0.05f64..1.0, snap_kind()), n))
    });
    verts
        .prop_map(move |((cx, cy), vs)| {
            let cx = e.xmin + cx * w;
            let cy = e.ymin + cy * h;
            let mut pts: Vec<(f64, f64, Snap)> = vs
                .iter()
                .enumerate()
                .map(|(i, &(jit, rad, s))| {
                    let ang = (i as f64 + jit * 0.9) / vs.len() as f64 * std::f64::consts::TAU;
                    (cx + ang.cos() * rad * w * 0.6, cy + ang.sin() * rad * h * 0.6, s)
                })
                .collect();
            pts.sort_by(|a, b| {
                let aa = (a.1 - cy).atan2(a.0 - cx);
                let ab = (b.1 - cy).atan2(b.0 - cx);
                aa.partial_cmp(&ab).unwrap()
            });
            let mut ring: Vec<(f64, f64)> = pts.iter().map(|&(x, y, s)| apply_snap(&g, x, y, s)).collect();
            ring.dedup();
            if ring.first() != ring.last() {
                let f = ring[0];
                ring.push(f);
            }
            (g.clone(), ring)
        })
        .prop_filter("simple ring with area", |(_, ring)| {
            ring.len() >= 4 && shoelace(ring).abs() > 1e-12 && is_simple(ring)
        })
}

fn ring() -> impl Strategy<Value = (G, Vec<(f64, f64)>)> {
    grid().prop_flat_map(ring_on)
}

// ------------------------------------------------------------- geometry --

fn shoelace(r: &[(f64, f64)]) -> f64 {
    let mut a = 0.0;
    for i in 0..r.len() - 1 {
        a += r[i].0 * r[i + 1].1 - r[i + 1].0 * r[i].1;
    }
    a / 2.0
}

fn orient(p: (f64, f64), q: (f64, f64), r: (f64, f64)) -> f64 {
    (q.0 - p.0) * (r.1 - p.1) - (q.1 - p.1) * (r.0 - p.0)
}

fn on_segment(p: (f64, f64), q: (f64, f64), r: (f64, f64)) -> bool {
    r.0 >= p.0.min(q.0) && r.0 <= p.0.max(q.0) && r.1 >= p.1.min(q.1) && r.1 <= p.1.max(q.1)
}

/// Segments ab and cd intersect or touch anywhere (closed test).
fn seg_touch(a: (f64, f64), b: (f64, f64), c: (f64, f64), d: (f64, f64)) -> bool {
    let o1 = orient(a, b, c);
    let o2 = orient(a, b, d);
    let o3 = orient(c, d, a);
    let o4 = orient(c, d, b);
    if (o1 > 0.0) != (o2 > 0.0) && (o3 > 0.0) != (o4 > 0.0) && o1 != 0.0 && o2 != 0.0 && o3 != 0.0 && o4 != 0.0 {
        return true;
    }
    (o1 == 0.0 && on_segment(a, b, c))
        || (o2 == 0.0 && on_segment(a, b, d))
        || (o3 == 0.0 && on_segment(c, d, a))
        || (o4 == 0.0 && on_segment(c, d, b))
}

/// A strictly simple ring: no repeated vertex, no two non-adjacent edges
/// that intersect or touch, no adjacent edges that fold back on each
/// other. Lattice snapping readily makes self-touching rings (a vertex
/// landing on another edge, two lobes meeting at a corner); those are
/// invalid polygons and the crate gives no defined result for them, so
/// the generator filters them.
fn is_simple(r: &[(f64, f64)]) -> bool {
    let n = r.len() - 1;
    for i in 0..n {
        for j in i + 1..n {
            if r[i] == r[j] {
                return false;
            }
        }
    }
    for i in 0..n {
        for j in i + 1..n {
            let adjacent = j == i + 1 || (i == 0 && j == n - 1);
            if adjacent {
                // fold-back: collinear and overlapping beyond the shared vertex
                let (a, b) = (r[i], r[i + 1]);
                let (c, d) = (r[j], r[j + 1]);
                if orient(a, b, c) == 0.0 && orient(a, b, d) == 0.0 {
                    let shared = if j == i + 1 { b } else { a };
                    let other_ab = if shared == b { a } else { b };
                    let other_cd = if shared == c { d } else { c };
                    let dot = (other_ab.0 - shared.0) * (other_cd.0 - shared.0)
                        + (other_ab.1 - shared.1) * (other_cd.1 - shared.1);
                    if dot > 0.0 {
                        return false;
                    }
                }
            } else if seg_touch(r[i], r[i + 1], r[j], r[j + 1]) {
                return false;
            }
        }
    }
    true
}

/// Sutherland-Hodgman clip of a ring against a rectangle; the area of the
/// result is the area of polygon-intersect-rectangle for any simple ring.
fn clip_to(r: &[(f64, f64)], e: &Ext) -> Vec<(f64, f64)> {
    let mut out: Vec<(f64, f64)> = r[..r.len() - 1].to_vec();
    type Inside = fn(&(f64, f64), f64) -> bool;
    type Isect = fn((f64, f64), (f64, f64), f64) -> (f64, f64);
    let planes: [(Inside, Isect, f64); 4] = [
        (
            |p, v| p.0 >= v,
            |a, b, v| {
                let t = (v - a.0) / (b.0 - a.0);
                (v, a.1 + t * (b.1 - a.1))
            },
            e.xmin,
        ),
        (
            |p, v| p.0 <= v,
            |a, b, v| {
                let t = (v - a.0) / (b.0 - a.0);
                (v, a.1 + t * (b.1 - a.1))
            },
            e.xmax,
        ),
        (
            |p, v| p.1 >= v,
            |a, b, v| {
                let t = (v - a.1) / (b.1 - a.1);
                (a.0 + t * (b.0 - a.0), v)
            },
            e.ymin,
        ),
        (
            |p, v| p.1 <= v,
            |a, b, v| {
                let t = (v - a.1) / (b.1 - a.1);
                (a.0 + t * (b.0 - a.0), v)
            },
            e.ymax,
        ),
    ];
    for (inside, isect, v) in planes.iter() {
        if out.is_empty() {
            break;
        }
        let input = out.clone();
        out.clear();
        for i in 0..input.len() {
            let cur = input[i];
            let prev = input[(i + input.len() - 1) % input.len()];
            let ci = inside(&cur, *v);
            let pi = inside(&prev, *v);
            if ci {
                if !pi {
                    out.push(isect(prev, cur, *v));
                }
                out.push(cur);
            } else if pi {
                out.push(isect(prev, cur, *v));
            }
        }
    }
    if !out.is_empty() {
        let f = out[0];
        out.push(f);
    }
    out
}

fn clipped_area(r: &[(f64, f64)], g: &G) -> f64 {
    let c = clip_to(r, &g.ext());
    if c.len() < 4 {
        0.0
    } else {
        shoelace(&c).abs()
    }
}

// -------------------------------------------------------------- helpers --

fn geom(r: &[(f64, f64)]) -> Geometry {
    Geometry::Polygon(Polygon::new(vec![r.iter().map(|&p| Coord::from(p)).collect()]))
}

/// Dense per-cell cover from the tables, row-major, `nrow * ncol`.
fn cover(res: &BurnResult, g: &G) -> Vec<f64> {
    let (nc, nr) = (g.ncol(), g.nrow());
    let mut m = vec![0.0; nc * nr];
    for x in &res.runs {
        for c in x.col_start..x.col_end {
            m[x.row as usize * nc + c as usize] += 1.0;
        }
    }
    for e in &res.edges {
        m[e.row as usize * nc + e.col as usize] += e.fraction as f64;
    }
    m
}

fn burn1(r: &[(f64, f64)], g: &G, opts: BurnOptions) -> BurnResult {
    burn(&[geom(r)], &g.spec, opts).unwrap()
}

fn fmt_ring(r: &[(f64, f64)]) -> String {
    let body: Vec<String> = r.iter().map(|p| format!("{:?} {:?}", p.0, p.1)).collect();
    format!("POLYGON (({}))", body.join(", "))
}

fn fmt_grid(g: &G) -> String {
    let e = g.spec.extent();
    format!("GridSpec::new({:?}, {:?}, {:?}, {:?}, {}, {})", e.xmin, e.ymin, e.xmax, e.ymax, g.ncol(), g.nrow())
}

// ----------------------------------------------------------- properties --

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 1500,
        max_shrink_iters: 4000,
        failure_persistence: Some(Box::new(proptest::test_runner::FileFailurePersistence::WithSource("proptest-regressions"))),
        .. ProptestConfig::default()
    })]

    #[test]
    fn coverage_conserves_area((g, ring) in ring()) {
        let r = burn1(&ring, &g, BurnOptions::default());
        let got = area(&r, &g.spec);
        let truth = clipped_area(&ring, &g);
        let cell = g.res() * g.res();
        // f32 fractions: absolute slack of a few ulp per edge cell
        let tol = 1e-6 * truth.max(cell) + r.edges.len() as f64 * cell * 2e-7;
        prop_assert!(r.notes.is_empty(), "notes: {:?}\n{}\n{}", r.notes, fmt_grid(&g), fmt_ring(&ring));
        prop_assert!((got - truth).abs() <= tol,
            "area {got} vs {truth} (diff {:+e}, tol {tol:e})\n{}\n{}", got - truth, fmt_grid(&g), fmt_ring(&ring));
    }

    #[test]
    fn tables_are_well_formed((g, ring) in ring()) {
        let r = burn1(&ring, &g, BurnOptions::default());
        let (nc, nr) = (g.ncol() as i32, g.nrow() as i32);
        for x in &r.runs {
            prop_assert!(x.row >= 0 && x.row < nr && x.col_start >= 0 && x.col_start < x.col_end && x.col_end <= nc,
                "bad run {x:?}\n{}\n{}", fmt_grid(&g), fmt_ring(&ring));
        }
        for e in &r.edges {
            prop_assert!(e.row >= 0 && e.row < nr && e.col >= 0 && e.col < nc, "bad edge {e:?}");
            prop_assert!(e.fraction > 0.0 && e.fraction <= 1.0 + 1e-6, "fraction out of range {e:?}\n{}\n{}", fmt_grid(&g), fmt_ring(&ring));
        }
        // a cell is either in a run or an edge, never both
        let (nc, _) = (g.ncol(), g.nrow());
        let mut seen = std::collections::HashSet::new();
        for x in &r.runs { for c in x.col_start..x.col_end { prop_assert!(seen.insert(x.row as usize * nc + c as usize), "cell twice (run) {x:?}"); } }
        for e in &r.edges { prop_assert!(seen.insert(e.row as usize * nc + e.col as usize), "cell twice (edge) {e:?}\n{}\n{}", fmt_grid(&g), fmt_ring(&ring)); }
    }

    #[test]
    fn invariant_under_orientation_and_start((g, ring) in ring(), k in 0usize..8) {
        let base = cover(&burn1(&ring, &g, BurnOptions::default()), &g);
        let rev: Vec<_> = ring.iter().rev().cloned().collect();
        let got = cover(&burn1(&rev, &g, BurnOptions::default()), &g);
        let worst = base.iter().zip(&got).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        prop_assert!(worst <= 1e-5, "reversed ring differs by {worst}\n{}\n{}", fmt_grid(&g), fmt_ring(&ring));
        let n = ring.len() - 1;
        let s = k % n;
        let mut rot: Vec<_> = ring[s..n].to_vec();
        rot.extend_from_slice(&ring[..s]);
        rot.push(rot[0]);
        let got = cover(&burn1(&rot, &g, BurnOptions::default()), &g);
        let worst = base.iter().zip(&got).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        prop_assert!(worst <= 1e-5, "start-rotated ring differs by {worst}\n{}\n{}", fmt_grid(&g), fmt_ring(&ring));
    }

    #[test]
    fn matches_exact_per_cell_oracle((g, ring) in grid_up_to(16).prop_flat_map(ring_on)) {
        let m = cover(&burn1(&ring, &g, BurnOptions::default()), &g);
        let (nc, nr) = (g.ncol(), g.nrow());
        let e = g.ext();
        let res = g.res();
        for row in 0..nr { for col in 0..nc {
            let cell = Ext { xmin: e.xmin + col as f64 * res, xmax: e.xmin + (col + 1) as f64 * res,
                             ymax: e.ymax - row as f64 * res, ymin: e.ymax - (row + 1) as f64 * res };
            let c = clip_to(&ring, &cell);
            let truth = if c.len() < 4 { 0.0 } else { shoelace(&c).abs() / (res * res) };
            let got = m[row * nc + col];
            prop_assert!((got - truth).abs() <= 2e-6,
                "row {row} col {col}: burned {got} vs clipped {truth}\n{}\n{}", fmt_grid(&g), fmt_ring(&ring));
        } }
    }

    #[test]
    fn approx_agrees_with_coverage_on_full_and_empty_cells((g, ring) in ring()) {
        let cov = cover(&burn1(&ring, &g, BurnOptions::default()), &g);
        let apx = cover(&burn1(&ring, &g, BurnOptions::approx()), &g);
        for (i, (c, a)) in cov.iter().zip(&apx).enumerate() {
            if *c >= 1.0 - 1e-6 {
                prop_assert!(*a > 0.0, "cell {i} fully covered but absent in approx\n{}\n{}", fmt_grid(&g), fmt_ring(&ring));
            }
            if *c == 0.0 {
                prop_assert!(*a == 0.0, "cell {i} untouched in coverage but present in approx\n{}\n{}", fmt_grid(&g), fmt_ring(&ring));
            }
        }
    }

    #[test]
    fn doubling_resolution_aggregates_exactly((g, ring) in ring()) {
        let m1 = cover(&burn1(&ring, &g, BurnOptions::default()), &g);
        let e = g.spec.extent();
        let g2 = G { spec: GridSpec::new(e.xmin, e.ymin, e.xmax, e.ymax, 2 * g.ncol() as u32, 2 * g.nrow() as u32) };
        let m2 = cover(&burn1(&ring, &g2, BurnOptions::default()), &g2);
        let (nc, nr) = (g.ncol(), g.nrow());
        for row in 0..nr { for col in 0..nc {
            let agg = (m2[(2 * row) * 2 * nc + 2 * col] + m2[(2 * row) * 2 * nc + 2 * col + 1]
                + m2[(2 * row + 1) * 2 * nc + 2 * col] + m2[(2 * row + 1) * 2 * nc + 2 * col + 1]) / 4.0;
            let got = m1[row * nc + col];
            prop_assert!((got - agg).abs() <= 2e-6, "row {row} col {col}: 1x {got} vs 2x aggregate {agg}\n{}\n{}", fmt_grid(&g), fmt_ring(&ring));
        } }
    }
}
