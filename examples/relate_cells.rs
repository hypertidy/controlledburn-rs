//! Prototype: DE-9IM of a polygon against every grid cell, from a finished
//! burn plus the geometry, checked cell by cell against `geo::Relate`.
//!
//! The burn is not changed. Two levels:
//!
//! 1. From the tables alone: II (2 if the cell has a record), EI (F for a
//!    run cell), EE (2), and IE / BE (F only if the polygon lies inside
//!    the closed cell, a bounding-box test).
//! 2. With the geometry, only for candidate cells: BI, BB, IB and EB.
//!    A cell the boundary touches without a record (f = 0, contact on a
//!    wall or a corner) must neighbour a cell the boundary passes through,
//!    so the candidates come from the tables: edge cells, run cells with a
//!    neighbour that is not a run, and the 8 neighbours of both (plus the
//!    grid's border cells, whose outer neighbours have no records). Every
//!    other cell has no boundary contact, and its matrix is level 1 with
//!    BI = BB = F and IB / EB from its record (run: IB = 1, EB = F; none:
//!    IB = F, EB = 1).
//!
//! Level 2 per candidate: each nearby boundary segment is tested against
//! the cell box with exact orientation tests (separating axes); BB is 1 for
//! a segment lying along a wall, 0 for point contact. IB / EB classify the
//! arcs of the cell's outline between boundary contacts: each arc lies
//! wholly inside or outside the polygon, so one point-in-polygon test at
//! its middle decides it.
//!
//! `cargo run --release --example relate_cells`. See design/de9im-cells.md.

use std::collections::{HashMap, HashSet};
use std::time::Instant;

use controlledburn::{burn, BurnOptions, BurnResult, Coord, Geometry, GridSpec, Polygon};
use geo::coordinate_position::CoordPos;
use geo::dimensions::Dimensions;
use geo::kernels::{Kernel, Orientation, RobustKernel};
use geo::{CoordinatePosition, Relate, Validation};

type GC = geo::Coord<f64>;

#[derive(Clone, Copy, PartialEq, Debug)]
enum Kind {
    Run,
    Edge,
    None,
}

#[derive(Clone, Copy)]
struct CellBox {
    x0: f64,
    y0: f64,
    x1: f64,
    y1: f64,
}

impl CellBox {
    fn corners(&self) -> [GC; 4] {
        [
            GC { x: self.x0, y: self.y0 },
            GC { x: self.x1, y: self.y0 },
            GC { x: self.x1, y: self.y1 },
            GC { x: self.x0, y: self.y1 },
        ]
    }
    fn strictly_inside(&self, p: GC) -> bool {
        p.x > self.x0 && p.x < self.x1 && p.y > self.y0 && p.y < self.y1
    }
    fn polygon(&self) -> geo::Polygon<f64> {
        geo::Rect::new(GC { x: self.x0, y: self.y0 }, GC { x: self.x1, y: self.y1 }).to_polygon()
    }
}

/// Unit cells on [0, n] x [0, n]: every wall is an exact integer.
fn cell_box(n: u32, row: i32, col: i32) -> CellBox {
    let top = n as f64 - row as f64;
    CellBox { x0: col as f64, y0: top - 1.0, x1: col as f64 + 1.0, y1: top }
}

fn orient(p: GC, q: GC, r: GC) -> i8 {
    match RobustKernel::orient2d(p, q, r) {
        Orientation::CounterClockwise => 1,
        Orientation::Clockwise => -1,
        Orientation::Collinear => 0,
    }
}

/// Does segment pq meet the open box? Separating axes: x, y and the
/// segment's normal (corners on the line count as separated).
fn meets_open(p: GC, q: GC, b: &CellBox) -> bool {
    if p.x.max(q.x) <= b.x0 || p.x.min(q.x) >= b.x1 || p.y.max(q.y) <= b.y0 || p.y.min(q.y) >= b.y1 {
        return false;
    }
    let o = b.corners().map(|c| orient(p, q, c));
    !(o.iter().all(|&v| v >= 0) || o.iter().all(|&v| v <= 0))
}

/// Does segment pq meet the closed box?
fn meets_closed(p: GC, q: GC, b: &CellBox) -> bool {
    if p.x.max(q.x) < b.x0 || p.x.min(q.x) > b.x1 || p.y.max(q.y) < b.y0 || p.y.min(q.y) > b.y1 {
        return false;
    }
    let o = b.corners().map(|c| orient(p, q, c));
    !(o.iter().all(|&v| v > 0) || o.iter().all(|&v| v < 0))
}

/// The part of pq lying along a wall with positive length, as an interval
/// of the outline's arc parameter, if any.
fn along_wall(p: GC, q: GC, b: &CellBox) -> Option<(f64, f64)> {
    let (w, h) = (b.x1 - b.x0, b.y1 - b.y0);
    if p.y == q.y && (p.y == b.y0 || p.y == b.y1) {
        let (lo, hi) = (p.x.min(q.x).max(b.x0), p.x.max(q.x).min(b.x1));
        if lo < hi {
            return Some(if p.y == b.y0 { (lo - b.x0, hi - b.x0) } else { (w + h + (b.x1 - hi), w + h + (b.x1 - lo)) });
        }
    }
    if p.x == q.x && (p.x == b.x0 || p.x == b.x1) {
        let (lo, hi) = (p.y.min(q.y).max(b.y0), p.y.max(q.y).min(b.y1));
        if lo < hi {
            return Some(if p.x == b.x1 {
                (w + (lo - b.y0), w + (hi - b.y0))
            } else {
                (2.0 * w + h + (b.y1 - hi), 2.0 * w + h + (b.y1 - lo))
            });
        }
    }
    None
}

/// Points where pq (not along a wall) meets the outline, as arc parameters.
fn wall_points(p: GC, q: GC, b: &CellBox, out: &mut Vec<f64>) {
    let (w, h) = (b.x1 - b.x0, b.y1 - b.y0);
    let mut hit = |s: f64| out.push(s);
    // vertical walls
    for (x, base, up) in [(b.x1, w, true), (b.x0, 2.0 * w + h, false)] {
        if (p.x - x) * (q.x - x) <= 0.0 && p.x != q.x {
            let t = (x - p.x) / (q.x - p.x);
            let y = if t == 0.0 {
                p.y
            } else if t == 1.0 {
                q.y
            } else {
                p.y + t * (q.y - p.y)
            };
            if y >= b.y0 && y <= b.y1 {
                hit(if up { base + (y - b.y0) } else { base + (b.y1 - y) });
            }
        } else if p.x == q.x && p.x == x {
            // along the wall's line but outside its span: a touch at an end
            for c in [p, q] {
                if c.y >= b.y0 && c.y <= b.y1 {
                    hit(if up { base + (c.y - b.y0) } else { base + (b.y1 - c.y) });
                }
            }
        }
    }
    // horizontal walls
    for (y, base, right) in [(b.y0, 0.0, true), (b.y1, w + h, false)] {
        if (p.y - y) * (q.y - y) <= 0.0 && p.y != q.y {
            let t = (y - p.y) / (q.y - p.y);
            let x = if t == 0.0 {
                p.x
            } else if t == 1.0 {
                q.x
            } else {
                p.x + t * (q.x - p.x)
            };
            if x >= b.x0 && x <= b.x1 {
                hit(if right { base + (x - b.x0) } else { base + (b.x1 - x) });
            }
        } else if p.y == q.y && p.y == y {
            for c in [p, q] {
                if c.x >= b.x0 && c.x <= b.x1 {
                    hit(if right { base + (c.x - b.x0) } else { base + (b.x1 - c.x) });
                }
            }
        }
    }
}

/// The point at arc parameter s on the outline (counterclockwise from the
/// bottom-left corner).
fn at(b: &CellBox, s: f64) -> GC {
    let (w, h) = (b.x1 - b.x0, b.y1 - b.y0);
    if s < w {
        GC { x: b.x0 + s, y: b.y0 }
    } else if s < w + h {
        GC { x: b.x1, y: b.y0 + (s - w) }
    } else if s < 2.0 * w + h {
        GC { x: b.x1 - (s - w - h), y: b.y1 }
    } else {
        GC { x: b.x0, y: b.y1 - (s - 2.0 * w - h) }
    }
}

struct Prepared<'a> {
    segs: Vec<(GC, GC)>,
    bins: HashMap<(i32, i32), Vec<usize>>,
    bbox: [f64; 4],
    geo: &'a geo::Polygon<f64>,
}

fn prepare(gp: &geo::Polygon<f64>, n: u32) -> Prepared<'_> {
    let mut segs = Vec::new();
    for ring in std::iter::once(gp.exterior()).chain(gp.interiors()) {
        for l in ring.lines() {
            if l.start != l.end {
                segs.push((l.start, l.end));
            }
        }
    }
    let mut bins: HashMap<(i32, i32), Vec<usize>> = HashMap::new();
    let cell = |x: f64| x.floor() as i32;
    let row = |y: f64| n as i32 - 1 - y.floor() as i32;
    for (i, (p, q)) in segs.iter().enumerate() {
        let (c0, c1) = (cell(p.x.min(q.x)) - 1, cell(p.x.max(q.x)) + 1);
        let (r0, r1) = (row(p.y.max(q.y)) - 1, row(p.y.min(q.y)) + 1);
        for r in r0.max(0)..=r1.min(n as i32 - 1) {
            for c in c0.max(0)..=c1.min(n as i32 - 1) {
                bins.entry((r, c)).or_default().push(i);
            }
        }
    }
    let ext = gp.exterior();
    let bbox = ext.0.iter().fold([f64::MAX, f64::MAX, f64::MIN, f64::MIN], |a, c| {
        [a[0].min(c.x), a[1].min(c.y), a[2].max(c.x), a[3].max(c.y)]
    });
    Prepared { segs, bins, bbox, geo: gp }
}

fn kinds(r: &BurnResult) -> HashMap<(i32, i32), Kind> {
    let mut k = HashMap::new();
    for run in &r.runs {
        for c in run.col_start..run.col_end {
            k.insert((run.row, c), Kind::Run);
        }
    }
    for e in &r.edges {
        k.insert((e.row, e.col), Kind::Edge);
    }
    k
}

fn candidates(k: &HashMap<(i32, i32), Kind>, n: u32) -> HashSet<(i32, i32)> {
    let n = n as i32;
    let kind = |r: i32, c: i32| k.get(&(r, c)).copied().unwrap_or(Kind::None);
    let mut seeds = Vec::new();
    for (&(r, c), &kd) in k {
        let boundaryish = kd == Kind::Edge
            || (-1..=1).any(|dr| (-1..=1).any(|dc| (dr, dc) != (0, 0) && kind(r + dr, c + dc) != Kind::Run));
        if boundaryish {
            seeds.push((r, c));
        }
    }
    let mut out = HashSet::new();
    for (r, c) in seeds {
        for dr in -1..=1 {
            for dc in -1..=1 {
                let (rr, cc) = (r + dr, c + dc);
                if rr >= 0 && rr < n && cc >= 0 && cc < n {
                    out.insert((rr, cc));
                }
            }
        }
    }
    for i in 0..n {
        out.extend([(0, i), (n - 1, i), (i, 0), (i, n - 1)]);
    }
    out
}

/// Exact candidates from the geometry: every cell whose closed box a
/// boundary segment meets. No record is needed, so no tolerance enters.
fn candidates_exact(p: &Prepared, n: u32) -> HashSet<(i32, i32)> {
    p.bins
        .iter()
        .filter(|(&(r, c), ids)| ids.iter().any(|&i| meets_closed(p.segs[i].0, p.segs[i].1, &cell_box(n, r, c))))
        .map(|(&cell, _)| cell)
        .collect()
}

/// Level 1 entries that need no geometry beyond the polygon's bbox.
fn level1(kd: Kind, b: &CellBox, bbox: [f64; 4]) -> [u8; 9] {
    let inside = bbox[0] >= b.x0 && bbox[2] <= b.x1 && bbox[1] >= b.y0 && bbox[3] <= b.y1;
    let ii = if kd == Kind::None { b'F' } else { b'2' };
    let ei = if kd == Kind::Run { b'F' } else { b'2' };
    let (ie, be) = if inside { (b'F', b'F') } else { (b'2', b'1') };
    // II IB IE / BI BB BE / EI EB EE; IB, BI, BB, EB filled later
    [ii, b'?', ie, b'?', b'?', be, ei, b'?', b'2']
}

fn relate_cell(p: &Prepared, kd: Kind, b: &CellBox, cell: (i32, i32), candidate: bool, exact: bool) -> [u8; 9] {
    let mut m = level1(kd, b, p.bbox);
    if !candidate {
        m[1] = if kd == Kind::Run { b'1' } else { b'F' };
        m[3] = b'F';
        m[4] = b'F';
        m[7] = if kd == Kind::Run { b'F' } else { b'1' };
        return m;
    }
    let (mut bi, mut bb_line, mut bb_point) = (false, false, false);
    let mut spans: Vec<(f64, f64)> = Vec::new();
    let mut points: Vec<f64> = Vec::new();
    for &i in p.bins.get(&cell).map(|v| v.as_slice()).unwrap_or(&[]) {
        let (s, e) = p.segs[i];
        if !meets_closed(s, e, b) {
            continue;
        }
        bi |= meets_open(s, e, b);
        if let Some(iv) = along_wall(s, e, b) {
            bb_line = true;
            spans.push(iv);
        } else if !(b.strictly_inside(s) && b.strictly_inside(e)) {
            bb_point = true;
            wall_points(s, e, b, &mut points);
        }
    }
    m[3] = if bi { b'1' } else { b'F' };
    if exact {
        // The record says whether coverage is above the burn's tolerance
        // (1e-6), not whether it is above 0. A boundary through the open
        // cell leaves area on both sides; without one, the open cell is
        // wholly inside or outside, and its centre says which.
        let (ii, ei) = if bi {
            (true, true)
        } else {
            let c = GC { x: (b.x0 + b.x1) / 2.0, y: (b.y0 + b.y1) / 2.0 };
            let inside = p.geo.coordinate_position(&c) == CoordPos::Inside;
            (inside, !inside)
        };
        m[0] = if ii { b'2' } else { b'F' };
        m[6] = if ei { b'2' } else { b'F' };
    }
    m[4] = if bb_line {
        b'1'
    } else if bb_point {
        b'0'
    } else {
        b'F'
    };
    // IB / EB from the arcs of the outline between contacts.
    let perim = 2.0 * ((b.x1 - b.x0) + (b.y1 - b.y0));
    spans.extend(points.iter().map(|&s| (s % perim, s % perim)));
    let (mut ib, mut eb) = (false, false);
    let mut classify = |s: f64| match p.geo.coordinate_position(&at(b, s.rem_euclid(perim))) {
        CoordPos::Inside => ib = true,
        CoordPos::Outside => eb = true,
        CoordPos::OnBoundary => {} // an arc too short for f64; counted as a failure below
    };
    if spans.is_empty() {
        classify(0.0);
    } else {
        spans.sort_by(|a, c| a.0.partial_cmp(&c.0).unwrap());
        let mut end = spans[0].1;
        let first = spans[0].0;
        for &(a, z) in &spans[1..] {
            if a > end {
                classify((end + a) / 2.0);
            }
            end = end.max(z);
        }
        if perim + first > end {
            classify((end + perim + first) / 2.0);
        }
    }
    m[1] = if ib { b'1' } else { b'F' };
    m[7] = if eb { b'1' } else { b'F' };
    m
}

fn oracle(gp: &geo::Polygon<f64>, b: &CellBox) -> [u8; 9] {
    let im = gp.relate(&b.polygon());
    let pos = [CoordPos::Inside, CoordPos::OnBoundary, CoordPos::Outside];
    let mut m = [b'F'; 9];
    for (i, a) in pos.iter().enumerate() {
        for (j, c) in pos.iter().enumerate() {
            m[3 * i + j] = match im.get(*a, *c) {
                Dimensions::Empty => b'F',
                Dimensions::ZeroDimensional => b'0',
                Dimensions::OneDimensional => b'1',
                Dimensions::TwoDimensional => b'2',
            };
        }
    }
    m
}

// ------------------------------------------------------------------ cases

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn to_cb(gp: &geo::Polygon<f64>) -> Geometry {
    let ring = |r: &geo::LineString<f64>| r.0.iter().map(|c| Coord::new(c.x, c.y)).collect::<Vec<_>>();
    let mut rings = vec![ring(gp.exterior())];
    rings.extend(gp.interiors().iter().map(ring));
    Geometry::Polygon(Polygon::new(rings))
}

fn poly(ext: Vec<(f64, f64)>, holes: Vec<Vec<(f64, f64)>>) -> geo::Polygon<f64> {
    geo::Polygon::new(ext.into(), holes.into_iter().map(Into::into).collect())
}

/// A random star; `snap` rounds vertices to that step (0 = no rounding).
fn star(rng: &mut Lcg, n: u32, snap: f64) -> geo::Polygon<f64> {
    let span = n as f64;
    let (cx, cy) = (span * (0.15 + 0.7 * rng.next()), span * (0.15 + 0.7 * rng.next()));
    let r = span * (0.05 + 0.25 * rng.next());
    let k = 5 + (rng.next() * 25.0) as usize;
    let mut pts: Vec<(f64, f64)> = (0..k)
        .map(|i| {
            let a = std::f64::consts::TAU * (i as f64 + 0.5 * rng.next()) / k as f64;
            let rr = r * (0.3 + 0.7 * rng.next());
            let (x, y) = (cx + rr * a.cos(), cy + rr * a.sin());
            if snap > 0.0 {
                ((x / snap).round() * snap, (y / snap).round() * snap)
            } else {
                (x, y)
            }
        })
        .collect();
    pts.dedup();
    pts.push(pts[0]);
    // Half get a hole: a smaller star about the same centre, well inside
    // (the outer star holds the disc of radius 0.18 r).
    let mut holes = vec![];
    if rng.next() < 0.5 {
        let k = 4 + (rng.next() * 8.0) as usize;
        let mut h: Vec<(f64, f64)> = (0..k)
            .map(|i| {
                let a = std::f64::consts::TAU * (i as f64 + 0.5 * rng.next()) / k as f64;
                let rr = r * (0.05 + 0.1 * rng.next());
                let (x, y) = (cx + rr * a.cos(), cy + rr * a.sin());
                if snap > 0.0 {
                    ((x / snap).round() * snap, (y / snap).round() * snap)
                } else {
                    (x, y)
                }
            })
            .collect();
        h.dedup();
        if h.len() >= 3 {
            h.push(h[0]);
            holes.push(h);
        }
    }
    poly(pts, holes)
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<(f64, f64)> {
    vec![(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]
}

fn fixed_cases() -> Vec<(&'static str, geo::Polygon<f64>)> {
    vec![
        ("wall-aligned square", poly(rect(2.0, 2.0, 6.0, 6.0), vec![])),
        ("equal to a cell", poly(rect(3.0, 3.0, 4.0, 4.0), vec![])),
        ("hole equal to a cell", poly(rect(1.0, 1.0, 6.0, 6.0), vec![rect(3.0, 3.0, 4.0, 4.0)])),
        ("strictly inside one cell", poly(vec![(3.2, 3.2), (3.8, 3.3), (3.5, 3.8), (3.2, 3.2)], vec![])),
        ("in one cell, vertices on walls", poly(vec![(3.0, 3.5), (4.0, 3.2), (3.6, 4.0), (3.0, 3.5)], vec![])),
        ("corner-touching diamond", poly(vec![(3.0, 2.0), (4.0, 3.0), (3.0, 4.0), (2.0, 3.0), (3.0, 2.0)], vec![])),
        ("diagonal through corners", poly(vec![(1.0, 1.0), (6.0, 6.0), (1.0, 6.0), (1.0, 1.0)], vec![])),
        ("crossing the extent", poly(rect(-2.5, 3.5, 4.5, 12.0), vec![])),
        ("along the extent", poly(rect(0.0, 0.0, 3.0, 8.0), vec![])),
        (
            "L-shape",
            poly(vec![(1.0, 1.0), (5.0, 1.0), (5.0, 2.5), (2.5, 2.5), (2.5, 6.0), (1.0, 6.0), (1.0, 1.0)], vec![]),
        ),
        (
            "hole touching outer wall",
            poly(rect(1.0, 1.0, 6.0, 6.0), vec![vec![(1.0, 3.0), (3.0, 2.0), (3.0, 4.0), (1.0, 3.0)]]),
        ),
    ]
}

#[derive(Default)]
struct Tally {
    polygons: usize,
    cells: usize,
    candidates: usize,
    wrong_cells: usize,
    wrong_entry: [usize; 9],
    burn_s: f64,
    relate_s: f64,
    oracle_s: f64,
    examples: Vec<String>,
}

fn check(name: &str, gp: &geo::Polygon<f64>, n: u32, exact: bool, t: &mut Tally) {
    let grid = GridSpec::new(0.0, 0.0, n as f64, n as f64, n, n);
    let t0 = Instant::now();
    let r = burn(&[to_cb(gp)], &grid, BurnOptions::default()).unwrap();
    t.burn_s += t0.elapsed().as_secs_f64();

    let t1 = Instant::now();
    let p = prepare(gp, n);
    let k = kinds(&r);
    let cand = if exact { candidates_exact(&p, n) } else { candidates(&k, n) };
    // every cell within two of the polygon's bbox; the rest are disjoint
    let (c0, c1) = ((p.bbox[0].floor() as i32 - 2).max(0), (p.bbox[2].ceil() as i32 + 2).min(n as i32 - 1));
    let (r0, r1) = (
        (n as i32 - 1 - p.bbox[3].ceil() as i32 - 2).max(0),
        (n as i32 - 1 - p.bbox[1].floor() as i32 + 2).min(n as i32 - 1),
    );
    let mut ours = Vec::new();
    for row in r0..=r1 {
        for col in c0..=c1 {
            let b = cell_box(n, row, col);
            let kd = k.get(&(row, col)).copied().unwrap_or(Kind::None);
            ours.push(((row, col), relate_cell(&p, kd, &b, (row, col), cand.contains(&(row, col)), exact)));
        }
    }
    t.relate_s += t1.elapsed().as_secs_f64();

    let t2 = Instant::now();
    t.polygons += 1;
    t.candidates += cand.len();
    for ((row, col), m) in ours {
        let b = cell_box(n, row, col);
        let truth = oracle(gp, &b);
        t.cells += 1;
        if m != truth {
            t.wrong_cells += 1;
            for i in 0..9 {
                t.wrong_entry[i] += usize::from(m[i] != truth[i]);
            }
            if t.examples.len() < 5 {
                use geo::{Area, BooleanOps};
                let f = gp.intersection(&b.polygon()).unsigned_area() / ((b.x1 - b.x0) * (b.y1 - b.y0));
                t.examples.push(format!(
                    "{name}: cell ({row}, {col}) ours {} geo {}, exact f = {f:.3e}, 1 - f = {:.3e}",
                    String::from_utf8_lossy(&m),
                    String::from_utf8_lossy(&truth),
                    1.0 - f
                ));
            }
        }
    }
    t.oracle_s += t2.elapsed().as_secs_f64();
}

fn report(label: &str, t: &Tally) {
    println!(
        "{label:<34} {:>5} polygons {:>8} cells {:>3} wrong  burn {:.3} s  relate {:.3} s  geo {:.3} s",
        t.polygons, t.cells, t.wrong_cells, t.burn_s, t.relate_s, t.oracle_s
    );
    if t.wrong_cells > 0 {
        let names = ["II", "IB", "IE", "BI", "BB", "BE", "EI", "EB", "EE"];
        let by: Vec<String> =
            (0..9).filter(|&i| t.wrong_entry[i] > 0).map(|i| format!("{} {}", names[i], t.wrong_entry[i])).collect();
        println!("    wrong by entry: {}", by.join(", "));
        for e in &t.examples {
            println!("    {e}");
        }
    }
}

fn main() {
    for (exact, title) in [
        (false, "A. candidates and II / EI from the records"),
        (true, "B. candidates from the geometry; II / EI from it on candidate cells"),
    ] {
        println!("{title}");
        let mut fixed = Tally::default();
        for (name, gp) in fixed_cases() {
            check(name, &gp, 8, exact, &mut fixed);
        }
        report("fixed cases (8 x 8)", &fixed);

        let mut rng = Lcg(42);
        for (label, snap) in
            [("random stars, float vertices", 0.0), ("stars on a 0.5 lattice", 0.5), ("stars on the cell lattice", 1.0)]
        {
            for n in [16u32, 64, 256] {
                let mut t = Tally::default();
                let mut tried = 0;
                while t.polygons < 200 && tried < 2000 {
                    tried += 1;
                    let gp = star(&mut rng, n, snap);
                    if gp.exterior().0.len() < 4 || !gp.is_valid() {
                        continue;
                    }
                    check(label, &gp, n, exact, &mut t);
                }
                report(&format!("{label}, {n} x {n}"), &t);
            }
        }
        println!();
    }
}
