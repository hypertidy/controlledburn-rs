//! `relate_cells` against `geo::Relate(polygon, cell box)` on every cell
//! near each polygon: the cell's relation if it has one, else COVERED for a
//! run cell and DISJOINT for the rest.

use std::collections::HashMap;

use controlledburn::{
    burn, relate_cells, BurnMode, BurnOptions, Coord, Geometry, GridSpec, IntersectionMatrix, Polygon,
};
use geo::coordinate_position::CoordPos;
use geo::dimensions::Dimensions;
use geo::{Relate, Validation};

type GP = geo::Polygon<f64>;

/// The burn's cell box: the last column and row stretched to the extent.
fn cell_box(g: &GridSpec, row: i64, col: i64) -> GP {
    let (dx, dy) = (g.dx(), g.dy());
    let x0 = g.xmin + col as f64 * dx;
    let x1 = if col == g.ncol as i64 - 1 { g.xmax } else { g.xmin + (col + 1) as f64 * dx };
    let y1 = g.ymax - row as f64 * dy;
    let y0 = if row == g.nrow as i64 - 1 { g.ymin } else { g.ymax - (row + 1) as f64 * dy };
    geo::Rect::new(geo::coord! { x: x0, y: y0 }, geo::coord! { x: x1, y: y1 }).to_polygon()
}

fn oracle(p: &geo::MultiPolygon<f64>, b: &GP) -> String {
    let im = p.relate(b);
    let pos = [CoordPos::Inside, CoordPos::OnBoundary, CoordPos::Outside];
    let mut s = String::new();
    for a in pos {
        for c in pos {
            s.push(match im.get(a, c) {
                Dimensions::Empty => 'F',
                Dimensions::ZeroDimensional => '0',
                Dimensions::OneDimensional => '1',
                Dimensions::TwoDimensional => '2',
            });
        }
    }
    s
}

fn to_cb(mp: &geo::MultiPolygon<f64>) -> Geometry {
    let ring = |r: &geo::LineString<f64>| r.0.iter().map(|c| Coord::new(c.x, c.y)).collect::<Vec<_>>();
    let polys: Vec<Polygon> =
        mp.0.iter()
            .map(|p| {
                let mut rings = vec![ring(p.exterior())];
                rings.extend(p.interiors().iter().map(ring));
                Polygon::new(rings)
            })
            .collect();
    if polys.len() == 1 {
        Geometry::Polygon(polys.into_iter().next().unwrap())
    } else {
        Geometry::MultiPolygon(polys)
    }
}

/// Checks every cell within two of each geometry's bbox; returns
/// (cells checked, wrong cells with a description).
fn check(geoms: &[geo::MultiPolygon<f64>], grid: &GridSpec, opts: BurnOptions) -> (usize, Vec<String>) {
    let cb: Vec<Geometry> = geoms.iter().map(to_cb).collect();
    let r = burn(&cb, grid, opts).unwrap();
    let rel = relate_cells(&r, &cb, grid, opts).unwrap();
    assert!(rel.windows(2).all(|w| (w[0].id, w[0].row, w[0].col) < (w[1].id, w[1].row, w[1].col)), "sorted, unique");
    let by: HashMap<(i32, i32, i32), IntersectionMatrix> =
        rel.iter().map(|c| ((c.id, c.row, c.col), c.matrix)).collect();
    let covered = |id: i32, row: i32, col: i32| {
        r.runs.iter().any(|x| x.id == id && x.row == row && x.col_start <= col && col < x.col_end)
    };
    let (dx, dy) = (grid.dx(), grid.dy());
    let (mut n, mut wrong) = (0, Vec::new());
    for (k, mp) in geoms.iter().enumerate() {
        let id = k as i32 + opts.id_offset;
        let pts = mp.0.iter().flat_map(|p| p.exterior().0.iter());
        let (x0, y0, x1, y1) = pts.fold((f64::MAX, f64::MAX, f64::MIN, f64::MIN), |a, c| {
            (a.0.min(c.x), a.1.min(c.y), a.2.max(c.x), a.3.max(c.y))
        });
        let clamp = |v: i64, hi: u32| v.clamp(0, hi as i64 - 1);
        let (c0, c1) = (
            clamp(((x0 - grid.xmin) / dx).floor() as i64 - 2, grid.ncol),
            clamp(((x1 - grid.xmin) / dx).floor() as i64 + 2, grid.ncol),
        );
        let (r0, r1) = (
            clamp(((grid.ymax - y1) / dy).floor() as i64 - 2, grid.nrow),
            clamp(((grid.ymax - y0) / dy).floor() as i64 + 2, grid.nrow),
        );
        for row in r0..=r1 {
            for col in c0..=c1 {
                let (row, col) = (row as i32, col as i32);
                let ours = by.get(&(id, row, col)).copied().unwrap_or(if covered(id, row, col) {
                    IntersectionMatrix::COVERED
                } else {
                    IntersectionMatrix::DISJOINT
                });
                let truth = oracle(mp, &cell_box(grid, row as i64, col as i64));
                n += 1;
                if ours.as_str() != truth {
                    wrong.push(format!("geometry {k}, cell ({row}, {col}): ours {ours}, geo {truth}"));
                }
            }
        }
    }
    (n, wrong)
}

fn assert_exact(label: &str, geoms: &[geo::MultiPolygon<f64>], grid: &GridSpec, opts: BurnOptions) -> usize {
    let (n, wrong) = check(geoms, grid, opts);
    assert!(
        wrong.is_empty(),
        "{label}: {} of {n} cells wrong, e.g.\n{}",
        wrong.len(),
        wrong[..wrong.len().min(8)].join("\n")
    );
    n
}

fn poly(ext: Vec<(f64, f64)>, holes: Vec<Vec<(f64, f64)>>) -> GP {
    geo::Polygon::new(ext.into(), holes.into_iter().map(Into::into).collect())
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<(f64, f64)> {
    vec![(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]
}

fn one(p: GP) -> geo::MultiPolygon<f64> {
    geo::MultiPolygon(vec![p])
}

#[test]
fn fixed_cases_on_a_unit_grid() {
    let cases = vec![
        one(poly(rect(2.0, 2.0, 6.0, 6.0), vec![])),
        one(poly(rect(3.0, 3.0, 4.0, 4.0), vec![])),
        one(poly(rect(1.0, 1.0, 6.0, 6.0), vec![rect(3.0, 3.0, 4.0, 4.0)])),
        one(poly(vec![(3.2, 3.2), (3.8, 3.3), (3.5, 3.8), (3.2, 3.2)], vec![])),
        one(poly(vec![(3.0, 3.5), (4.0, 3.2), (3.6, 4.0), (3.0, 3.5)], vec![])),
        one(poly(vec![(3.0, 2.0), (4.0, 3.0), (3.0, 4.0), (2.0, 3.0), (3.0, 2.0)], vec![])),
        one(poly(vec![(1.0, 1.0), (6.0, 6.0), (1.0, 6.0), (1.0, 1.0)], vec![])),
        one(poly(rect(-2.5, 3.5, 4.5, 12.0), vec![])),
        one(poly(rect(0.0, 0.0, 3.0, 8.0), vec![])),
        one(poly(rect(-1.0, -1.0, 9.0, 9.0), vec![])), // covers the grid
        one(poly(vec![(1.0, 1.0), (5.0, 1.0), (5.0, 2.5), (2.5, 2.5), (2.5, 6.0), (1.0, 6.0), (1.0, 1.0)], vec![])),
        one(poly(rect(1.0, 1.0, 6.0, 6.0), vec![vec![(1.0, 3.0), (3.0, 2.0), (3.0, 4.0), (1.0, 3.0)]])),
        // a hole whose corner touches the shell at a cell corner
        one(poly(rect(1.0, 1.0, 6.0, 6.0), vec![vec![(1.0, 1.0), (3.0, 2.0), (2.0, 3.0), (1.0, 1.0)]])),
        // two squares meeting at a corner, and a hole equal to a cell beside a part
        geo::MultiPolygon(vec![poly(rect(1.0, 1.0, 3.0, 3.0), vec![]), poly(rect(3.0, 3.0, 5.0, 5.0), vec![])]),
        geo::MultiPolygon(vec![
            poly(rect(1.0, 1.0, 4.0, 4.0), vec![rect(2.0, 2.0, 3.0, 3.0)]),
            poly(rect(5.5, 1.5, 7.0, 7.5), vec![]),
        ]),
    ];
    for c in &cases {
        assert!(c.is_valid(), "test case is not valid: {c:?}");
    }
    let grid = GridSpec::new(0.0, 0.0, 8.0, 8.0, 8, 8);
    let n = assert_exact("fixed", &cases, &grid, BurnOptions::default());
    assert!(n > 500);
}

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// A random star in [x0, x0 + span]^2, half with a hole, vertices rounded
/// to multiples of `snap` (0: not rounded).
fn star(rng: &mut Lcg, x0: f64, y0: f64, span: f64, snap: f64) -> GP {
    let (cx, cy) = (x0 + span * (0.15 + 0.7 * rng.next()), y0 + span * (0.15 + 0.7 * rng.next()));
    let r = span * (0.05 + 0.25 * rng.next());
    let round = |v: f64| if snap > 0.0 { (v / snap).round() * snap } else { v };
    let ring = |rng: &mut Lcg, k: usize, r0: f64, r1: f64| {
        let mut pts: Vec<(f64, f64)> = (0..k)
            .map(|i| {
                let a = std::f64::consts::TAU * (i as f64 + 0.5 * rng.next()) / k as f64;
                let rr = r0 + (r1 - r0) * rng.next();
                (round(cx + rr * a.cos()), round(cy + rr * a.sin()))
            })
            .collect();
        pts.dedup();
        pts.push(pts[0]);
        pts
    };
    let k = 5 + (rng.next() * 25.0) as usize;
    let ext = ring(rng, k, 0.3 * r, r);
    let mut holes = vec![];
    if rng.next() < 0.5 {
        let kh = 4 + (rng.next() * 8.0) as usize;
        let h = ring(rng, kh, 0.05 * r, 0.15 * r);
        if h.len() >= 4 {
            holes.push(h);
        }
    }
    poly(ext, holes)
}

fn stars(seed: u64, count: usize, grid: &GridSpec, snap: f64) -> Vec<geo::MultiPolygon<f64>> {
    let mut rng = Lcg(seed);
    let mut out = Vec::new();
    let (w, h) = (grid.xmax - grid.xmin, grid.ymax - grid.ymin);
    while out.len() < count {
        // a little past the extent too
        let p = star(&mut rng, grid.xmin - 0.1 * w, grid.ymin - 0.1 * h, 1.2 * w.min(h), snap);
        if p.exterior().0.len() >= 4 && p.is_valid() {
            out.push(one(p));
        }
    }
    out
}

#[test]
fn random_polygons_float_vertices() {
    for n in [16u32, 48] {
        let grid = GridSpec::new(0.0, 0.0, n as f64, n as f64, n, n);
        assert_exact("float", &stars(1 + n as u64, 60, &grid, 0.0), &grid, BurnOptions::default());
    }
}

#[test]
fn random_polygons_vertices_on_walls_and_corners() {
    for snap in [0.5, 1.0] {
        for n in [16u32, 48] {
            let grid = GridSpec::new(0.0, 0.0, n as f64, n as f64, n, n);
            assert_exact("lattice", &stars(7 + n as u64, 60, &grid, snap), &grid, BurnOptions::default());
        }
    }
}

#[test]
fn non_dyadic_grid() {
    let grid = GridSpec::new(-3.3, 1.7, 9.1, 12.4, 37, 29);
    assert_exact("non-dyadic", &stars(99, 60, &grid, 0.0), &grid, BurnOptions::default());
}

#[test]
fn multipolygons() {
    let grid = GridSpec::new(0.0, 0.0, 32.0, 32.0, 32, 32);
    let mut rng = Lcg(5);
    let mut mps = Vec::new();
    while mps.len() < 30 {
        let a = star(&mut rng, 0.0, 0.0, 15.0, 0.0);
        let b = star(&mut rng, 16.5, 16.5, 15.0, 0.5);
        let mp = geo::MultiPolygon(vec![a, b]);
        if mp.is_valid() {
            mps.push(mp);
        }
    }
    assert_exact("multipolygons", &mps, &grid, BurnOptions::default());
}

#[test]
fn any_burn_options() {
    let grid = GridSpec::new(0.0, 0.0, 24.0, 24.0, 24, 24);
    let geoms = stars(3, 30, &grid, 0.5);
    for opts in [
        BurnOptions::approx(),
        BurnOptions::parity(BurnMode::Coverage),
        BurnOptions::parity(BurnMode::Approx),
        BurnOptions::default().with_coalesce_runs(false),
        BurnOptions::default().with_id_offset(5),
    ] {
        assert_exact(&format!("{opts:?}"), &geoms, &grid, opts);
    }
}

#[test]
fn lines_points_and_empty_give_nothing() {
    let grid = GridSpec::new(0.0, 0.0, 4.0, 4.0, 4, 4);
    let geoms = vec![
        Geometry::LineString(vec![Coord::new(0.5, 0.5), Coord::new(3.5, 3.5)]),
        Geometry::Point(Coord::new(1.5, 1.5)),
        Geometry::MultiPolygon(vec![]),
        Geometry::Polygon(Polygon::new(vec![vec![
            Coord::new(0.5, 0.5),
            Coord::new(f64::NAN, 0.5),
            Coord::new(1.5, 1.5),
            Coord::new(0.5, 0.5),
        ]])),
    ];
    let r = burn(&geoms, &grid, BurnOptions::default()).unwrap();
    assert!(relate_cells(&r, &geoms, &grid, BurnOptions::default()).unwrap().is_empty());
}

#[test]
fn errors_as_burn() {
    let grid = GridSpec::new(0.0, 0.0, 4.0, 4.0, 4, 4);
    let r = burn(&[], &grid, BurnOptions::default()).unwrap();
    assert!(relate_cells(&r, &[], &GridSpec::new(0.0, 0.0, 0.0, 4.0, 4, 4), BurnOptions::default()).is_err());
    assert!(relate_cells(&r, &[], &grid, BurnOptions::default().with_id_offset(-1)).is_err());
}
