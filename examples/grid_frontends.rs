//! Coverage fractions on grids that are not a regular lattice, using the
//! burn engine unchanged: a front end moves each polygon into index space
//! (fractional column, row), the engine burns it on a lattice of unit
//! cells, and the cells come back as (row, col) of the original grid.
//!
//! Three grids, three front ends:
//!
//! - rotated affine: the inverse geotransform, vertices only (exact);
//! - rectilinear (1D coordinates with uneven spacing): each edge cut
//!   where it crosses a cell wall, then mapped (exact);
//! - curvilinear (a polar sector, 2D coordinates): each edge cut into
//!   short pieces, then mapped (approximate: see the note it prints).
//!
//! Each is checked against a world-space reference: every cell outline
//! clipped by the polygon in the polygon's own CRS. Run with
//!
//! ```text
//! cargo run --release --example grid_frontends
//! ```

use std::collections::BTreeMap;

use controlledburn::{burn, BurnOptions, Coord, Geometry, GridSpec, Polygon};

type Ring = Vec<[f64; 2]>;

/// A named grid, its polygons, and points per cell side for its outlines.
type Case<'a> = (&'a str, &'a dyn Grid, Vec<Vec<Ring>>, usize);

/// A grid that is not (necessarily) a regular lattice in its CRS.
trait Grid {
    /// `(ncol, nrow)`.
    fn shape(&self) -> (usize, usize);
    /// Position `(x, y)` to fractional `(col, row)`, row 0 at the top.
    fn to_index(&self, p: [f64; 2]) -> [f64; 2];
    /// Fractional `(col, row)` to position.
    fn to_world(&self, i: [f64; 2]) -> [f64; 2];
    /// Cut a ring's edges so that each piece maps to a straight segment in
    /// index space (exactly, or as nearly as the grid allows), then map.
    fn ring_to_index(&self, r: &Ring) -> Ring;
}

/// GDAL-order geotransform `[x0, dx, rx, y0, ry, dy]`, rotation allowed.
struct Affine {
    gt: [f64; 6],
    inv: [f64; 6],
    n: (usize, usize),
}

impl Affine {
    fn new(gt: [f64; 6], ncol: usize, nrow: usize) -> Self {
        let det = gt[1] * gt[5] - gt[2] * gt[4];
        let (a, b, c, d) = (gt[5] / det, -gt[2] / det, -gt[4] / det, gt[1] / det);
        let inv = [-a * gt[0] - b * gt[3], a, b, -c * gt[0] - d * gt[3], c, d];
        Affine { gt, inv, n: (ncol, nrow) }
    }
}

impl Grid for Affine {
    fn shape(&self) -> (usize, usize) {
        self.n
    }
    fn to_index(&self, p: [f64; 2]) -> [f64; 2] {
        let g = &self.inv;
        [g[0] + p[0] * g[1] + p[1] * g[2], g[3] + p[0] * g[4] + p[1] * g[5]]
    }
    fn to_world(&self, i: [f64; 2]) -> [f64; 2] {
        let g = &self.gt;
        [g[0] + i[0] * g[1] + i[1] * g[2], g[3] + i[0] * g[4] + i[1] * g[5]]
    }
    /// An affine map keeps straight edges straight: vertices are enough.
    fn ring_to_index(&self, r: &Ring) -> Ring {
        r.iter().map(|p| self.to_index(*p)).collect()
    }
}

/// Cell walls along x (ascending) and y (descending: row 0 at the top).
struct Rectilinear {
    xw: Vec<f64>,
    yw: Vec<f64>,
}

/// Fractional index of `v` among walls `w` (ascending or descending),
/// linear within a cell and extended past the ends by the end cells.
fn wall_index(w: &[f64], v: f64) -> f64 {
    let n = w.len() - 1;
    let asc = w[n] > w[0];
    // first wall beyond v, in the walls' direction
    let k = w.partition_point(|&b| if asc { b <= v } else { b >= v });
    let i = k.saturating_sub(1).min(n - 1);
    i as f64 + (v - w[i]) / (w[i + 1] - w[i])
}

fn wall_value(w: &[f64], i: f64) -> f64 {
    let n = w.len() - 1;
    let k = (i.floor().max(0.0) as usize).min(n - 1);
    w[k] + (i - k as f64) * (w[k + 1] - w[k])
}

impl Grid for Rectilinear {
    fn shape(&self) -> (usize, usize) {
        (self.xw.len() - 1, self.yw.len() - 1)
    }
    fn to_index(&self, p: [f64; 2]) -> [f64; 2] {
        [wall_index(&self.xw, p[0]), wall_index(&self.yw, p[1])]
    }
    fn to_world(&self, i: [f64; 2]) -> [f64; 2] {
        [wall_value(&self.xw, i[0]), wall_value(&self.yw, i[1])]
    }
    /// The map is affine within each cell, so cutting each edge where it
    /// crosses a wall makes the index-space ring exact.
    fn ring_to_index(&self, r: &Ring) -> Ring {
        let mut out = Vec::new();
        for s in r.windows(2) {
            let (a, b) = (s[0], s[1]);
            let mut ts = vec![0.0];
            for (axis, walls) in [(0, &self.xw), (1, &self.yw)] {
                let (lo, hi) = (a[axis].min(b[axis]), a[axis].max(b[axis]));
                for &w in walls.iter() {
                    if w > lo && w < hi {
                        ts.push((w - a[axis]) / (b[axis] - a[axis]));
                    }
                }
            }
            ts.sort_by(f64::total_cmp);
            for t in ts {
                out.push(self.to_index([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]));
            }
        }
        out.push(self.to_index(*r.last().unwrap()));
        out
    }
}

/// A polar sector: radius `r0..r1` across columns, angle `a0..a1` degrees
/// down rows. Its cells are annular sectors; the inverse is analytic here,
/// where a general curvilinear grid would search its coordinate arrays.
struct Sector {
    r0: f64,
    r1: f64,
    a0: f64,
    a1: f64,
    n: (usize, usize),
    /// Longest piece, in cells, an edge is cut into.
    step: f64,
}

impl Grid for Sector {
    fn shape(&self) -> (usize, usize) {
        self.n
    }
    fn to_index(&self, p: [f64; 2]) -> [f64; 2] {
        let r = p[0].hypot(p[1]);
        let a = p[1].atan2(p[0]).to_degrees();
        [(r - self.r0) / (self.r1 - self.r0) * self.n.0 as f64, (a - self.a0) / (self.a1 - self.a0) * self.n.1 as f64]
    }
    fn to_world(&self, i: [f64; 2]) -> [f64; 2] {
        let r = self.r0 + (self.r1 - self.r0) * i[0] / self.n.0 as f64;
        let a = (self.a0 + (self.a1 - self.a0) * i[1] / self.n.1 as f64).to_radians();
        [r * a.cos(), r * a.sin()]
    }
    /// Straight edges become curves; cut them into short pieces.
    fn ring_to_index(&self, r: &Ring) -> Ring {
        let mut out = Vec::new();
        for s in r.windows(2) {
            let (a, b) = (s[0], s[1]);
            let (ia, ib) = (self.to_index(a), self.to_index(b));
            let len = (ib[0] - ia[0]).hypot(ib[1] - ia[1]);
            let k = (len / self.step).ceil().max(1.0) as usize;
            for j in 0..k {
                let t = j as f64 / k as f64;
                out.push(self.to_index([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]));
            }
        }
        out.push(self.to_index(*r.last().unwrap()));
        out
    }
}

/// The front end: rings to index space (y negated, so rows run down),
/// burned on unit cells; fractions summed per cell with holes subtracted.
/// With `cut` false, vertices are mapped one by one (edges left straight in
/// index space), for comparison.
fn burn_index(g: &dyn Grid, polygon: &[Ring], cut: bool) -> BTreeMap<(usize, usize), f64> {
    let (w, h) = g.shape();
    let grid = GridSpec::new(0.0, -(h as f64), w as f64, 0.0, w as u32, h as u32);
    let map = |r: &Ring| if cut { g.ring_to_index(r) } else { r.iter().map(|p| g.to_index(*p)).collect() };
    let rings: Vec<Vec<Coord>> =
        polygon.iter().map(|r| map(r).iter().map(|p| Coord::new(p[0], -p[1])).collect()).collect();
    let geoms = vec![Geometry::Polygon(Polygon::new(rings))];
    let r = burn(&geoms, &grid, BurnOptions::coverage()).unwrap();
    let mut out = BTreeMap::new();
    for run in &r.runs {
        for c in run.col_start..run.col_end {
            out.insert((run.row as usize, c as usize), 1.0);
        }
    }
    for e in &r.edges {
        out.insert((e.row as usize, e.col as usize), e.fraction as f64);
    }
    out
}

/// The reference: each cell outline (`k` points a side) clipped by the
/// polygon in world space, as a fraction of the outline's area.
fn clip_world(g: &dyn Grid, polygon: &[Ring], k: usize) -> (BTreeMap<(usize, usize), f64>, f64) {
    let (w, h) = g.shape();
    let mut out = BTreeMap::new();
    let mut covered = 0.0;
    for row in 0..h {
        for col in 0..w {
            let (c, r) = (col as f64, row as f64);
            let mut q = Vec::new();
            for side in 0..4 {
                for j in 0..k {
                    let t = j as f64 / k as f64;
                    let i = match side {
                        0 => [c + t, r],
                        1 => [c + 1.0, r + t],
                        2 => [c + 1.0 - t, r + 1.0],
                        _ => [c, r + 1.0 - t],
                    };
                    q.push(g.to_world(i));
                }
            }
            let a = overlap(polygon, &q);
            if a > 0.0 {
                let cell = area(&q).abs();
                out.insert((row, col), a / cell);
                covered += a;
            }
        }
    }
    (out, covered)
}

fn area(p: &[[f64; 2]]) -> f64 {
    let n = p.len();
    (0..n).map(|i| p[i][0] * p[(i + 1) % n][1] - p[(i + 1) % n][0] * p[i][1]).sum::<f64>() / 2.0
}

/// Area of `polygon` inside the simple outline `q`: `q` fanned into
/// signed triangles, the polygon clipped by each (Sutherland-Hodgman).
fn overlap(polygon: &[Ring], q: &[[f64; 2]]) -> f64 {
    let mut total = 0.0;
    for i in 1..q.len() - 1 {
        let t = [q[0], q[i], q[i + 1]];
        let s = area(&t);
        if s == 0.0 {
            continue;
        }
        let t = if s > 0.0 { t } else { [q[0], q[i + 1], q[i]] };
        let mut a = 0.0;
        for (k, ring) in polygon.iter().enumerate() {
            let c = area(&clip(ring, &t)).abs();
            a += if k == 0 { c } else { -c };
        }
        total += s.signum() * a;
    }
    total.abs()
}

fn clip(subject: &Ring, t: &[[f64; 2]; 3]) -> Ring {
    let mut out = subject[..subject.len() - 1].to_vec();
    for i in 0..3 {
        let (a, b) = (t[i], t[(i + 1) % 3]);
        let side = |p: [f64; 2]| (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]);
        let input = std::mem::take(&mut out);
        for j in 0..input.len() {
            let (p, q) = (input[j], input[(j + 1) % input.len()]);
            let (sp, sq) = (side(p), side(q));
            if sp >= 0.0 {
                out.push(p);
            }
            if (sp >= 0.0) != (sq >= 0.0) {
                let u = sp / (sp - sq);
                out.push([p[0] + u * (q[0] - p[0]), p[1] + u * (q[1] - p[1])]);
            }
        }
        if out.is_empty() {
            break;
        }
    }
    out
}

/// A five-pointed star with a square hole.
fn star(c: [f64; 2], r: f64) -> Vec<Ring> {
    let mut ext: Ring = (0..10)
        .map(|k| {
            let a = (90.0 + 36.0 * k as f64).to_radians();
            let rr = if k % 2 == 0 { r } else { 0.45 * r };
            [c[0] + rr * a.cos(), c[1] + rr * a.sin()]
        })
        .collect();
    ext.push(ext[0]);
    let h = 0.15 * r;
    let hole = vec![
        [c[0] - h, c[1] - h],
        [c[0] - h, c[1] + h],
        [c[0] + h, c[1] + h],
        [c[0] + h, c[1] - h],
        [c[0] - h, c[1] - h],
    ];
    vec![ext, hole]
}

/// A long thin quadrilateral from `a` to `b`, `w` wide.
fn sliver(a: [f64; 2], b: [f64; 2], w: f64) -> Vec<Ring> {
    let l = (b[0] - a[0]).hypot(b[1] - a[1]);
    let (nx, ny) = (-(b[1] - a[1]) / l * w / 2.0, (b[0] - a[0]) / l * w / 2.0);
    vec![vec![
        [a[0] + nx, a[1] + ny],
        [a[0] - nx, a[1] - ny],
        [b[0] - nx, b[1] - ny],
        [b[0] + nx, b[1] + ny],
        [a[0] + nx, a[1] + ny],
    ]]
}

fn polygon_area(p: &[Ring]) -> f64 {
    area(&p[0][..p[0].len() - 1]).abs() - p[1..].iter().map(|r| area(&r[..r.len() - 1]).abs()).sum::<f64>()
}

/// Walls of `n` cells from `w0`, each `r` times wider than the last.
fn stretched(n: usize, w0: f64, d0: f64, r: f64) -> Vec<f64> {
    let mut w = vec![w0];
    let mut d = d0;
    for _ in 0..n {
        w.push(w.last().unwrap() + d);
        d *= r;
    }
    w
}

fn main() {
    let t = 30f64.to_radians();
    let rotated =
        Affine::new([0.0, 1000.0 * t.cos(), 1000.0 * t.sin(), 30_000.0, 1000.0 * t.sin(), -1000.0 * t.cos()], 40, 30);
    let mut yw = stretched(28, 0.0, 1500.0, 0.93);
    yw.reverse();
    let rect = Rectilinear { xw: stretched(36, 0.0, 300.0, 1.08), yw };
    let sector = Sector { r0: 8_000.0, r1: 40_000.0, a0: 100.0, a1: 10.0, n: (24, 30), step: 0.05 };

    let cases: Vec<Case> = vec![
        (
            "rotated affine",
            &rotated,
            vec![star([24_820.0, 27_010.0], 9_000.0), sliver([12_000.0, 18_000.0], [33_000.0, 34_000.0], 1_500.0)],
            1,
        ),
        (
            "rectilinear",
            &rect,
            vec![star([12_000.0, 9_000.0], 6_500.0), sliver([3_000.0, 2_000.0], [36_000.0, 14_000.0], 1_200.0)],
            1,
        ),
        (
            "polar sector",
            &sector,
            vec![star([8_200.0, 22_550.0], 9_000.0), sliver([-1_000.0, 12_000.0], [36_000.0, 10_000.0], 1_500.0)],
            64,
        ),
    ];

    println!(
        "{:16} {:>7} {:>16} {:>16} {:>18} {:>18}",
        "grid", "polygon", "area error, ref", "area error, burn", "largest cell diff", "vertices only"
    );
    for (name, g, polys, k) in cases {
        for (i, p) in polys.iter().enumerate() {
            let truth = polygon_area(p);
            let (reference, covered) = clip_world(g, p, k);
            let burned = burn_index(g, p, true);
            let verts = burn_index(g, p, false);
            let mut worst_verts: f64 = 0.0;
            for c in reference.keys().chain(verts.keys()) {
                let (a, b) = (reference.get(c).copied().unwrap_or(0.0), verts.get(c).copied().unwrap_or(0.0));
                worst_verts = worst_verts.max((a - b).abs());
            }
            // burned fractions weighted by each cell's world area
            let mut burned_area = 0.0;
            let mut worst: f64 = 0.0;
            let mut cells: BTreeMap<(usize, usize), (f64, f64)> = BTreeMap::new();
            for (c, f) in &reference {
                cells.entry(*c).or_default().0 = *f;
            }
            for (c, f) in &burned {
                cells.entry(*c).or_default().1 = *f;
            }
            for ((row, col), (fr, fb)) in cells {
                let (cc, rr) = (col as f64, row as f64);
                let q: Vec<[f64; 2]> = (0..4 * k)
                    .map(|j| {
                        let (side, t) = (j / k, (j % k) as f64 / k as f64);
                        g.to_world(match side {
                            0 => [cc + t, rr],
                            1 => [cc + 1.0, rr + t],
                            2 => [cc + 1.0 - t, rr + 1.0],
                            _ => [cc, rr + 1.0 - t],
                        })
                    })
                    .collect();
                burned_area += fb * area(&q).abs();
                worst = worst.max((fr - fb).abs());
            }
            println!(
                "{:16} {:>7} {:>16.1e} {:>16.1e} {:>18.4} {:>18.4}",
                name,
                i,
                covered / truth - 1.0,
                burned_area / truth - 1.0,
                worst,
                worst_verts
            );
        }
    }
    println!(
        "\nRotated and rectilinear: the burn matches the world-space reference to single precision.\n\
         Polar sector: the edges are followed closely, but a fraction of a unit cell is not a\n\
         fraction of the cell's area (the map is not affine within a cell), so cells differ by\n\
         about a percent while the total area stays close.\n\
         Vertices only (edges straight in index space) is wrong by whole cells on both\n\
         non-affine grids."
    );
}
