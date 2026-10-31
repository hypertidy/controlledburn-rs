//! Canonical tables for the C++ core.
//!
//! This crate is the canon for the burn algorithm; the C++ core
//! (hypertidy/controlledburn, `cpp/`) follows it. This test generates a
//! deterministic set of lattice-biased cases (rings with vertices on cell
//! corners, walls and centre lines, self-intersections from snapping,
//! holes, lines and points on walls) on non-dyadic grids, burns each in
//! parity mode (no run coalescing, as the C++ core emits) in both modes,
//! and compares with `fixtures/canon-tables.txt` record by record, with
//! fractions and lengths compared as bits.
//!
//! The same file is copied into the C++ repository as
//! `fixtures/canon-tables.txt`, where `cpp/tests/test_canon_tables.cpp`
//! checks the C++ core against it. After a deliberate change to the
//! algorithm, regenerate it with
//!
//!   CANON_WRITE=1 cargo test --release --test canon
//!
//! and port the change to the C++ core until its test passes on the new
//! file.
//!
//! Format, one record per line: `case <name>`, `grid <xmin> <ymin> <xmax>
//! <ymax> <ncol> <nrow>` (bounds as IEEE-754 bit patterns in hex), one
//! `wkb <hex>` per geometry (geometry k has id k), then for each of
//! `mode coverage` and `mode approx`: `run <row> <col_start> <col_end>
//! <id>`, `edge <row> <col> <f32 bits> <id>`, `line <row> <col> <f32 bits>
//! <id>` and `point <row> <col> <id>`, all 0-based with half-open runs;
//! `end` closes a case.

use controlledburn::*;
use std::fmt::Write as _;

struct Lcg(u64);
impl Lcg {
    fn u(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn pick<T: Copy>(&mut self, v: &[T]) -> T {
        v[((self.u() * v.len() as f64) as usize).min(v.len() - 1)]
    }
}

#[derive(Clone, Copy)]
struct G {
    x0: f64,
    y1: f64,
    r: f64,
    nc: u32,
    nr: u32,
}

impl G {
    fn spec(&self) -> GridSpec {
        GridSpec::new(
            self.x0,
            self.y1 - self.nr as f64 * self.r,
            self.x0 + self.nc as f64 * self.r,
            self.y1,
            self.nc,
            self.nr,
        )
    }
    /// A point inside (or a little outside) the grid, steered onto a seam:
    /// lattice positions are recomputed from the origin as a user would
    /// type them, so they disagree with the grid's own by an ulp.
    fn point(&self, g: &mut Lcg, cx: f64, cy: f64, rad: f64, ang: f64) -> Coord {
        let w = self.nc as f64 * self.r;
        let h = self.nr as f64 * self.r;
        let (x, y) = (cx + ang.cos() * rad * w * 0.6, cy + ang.sin() * rad * h * 0.6);
        let kx = ((x - self.x0) / self.r).round();
        let ky = ((self.y1 - y) / self.r).round();
        let (xl, yl) = (self.x0 + kx * self.r, self.y1 - ky * self.r);
        let (xc, yc) = (self.x0 + (kx + 0.5) * self.r, self.y1 - (ky + 0.5) * self.r);
        let c = match g.pick(&[0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6]) {
            0 => (x, y),
            1 => (xl, yl),
            2 => (xl, y),
            3 => (x, yl),
            4 => (x, yc),
            5 => (xc, y),
            _ => (xc, yc),
        };
        Coord::new(c.0, c.1)
    }
    fn ring(&self, g: &mut Lcg, cx: f64, cy: f64, scale: f64) -> Vec<Coord> {
        let n = 3 + (g.u() * 9.0) as usize;
        let mut v = Vec::with_capacity(n + 1);
        for i in 0..n {
            let ang = (i as f64 + g.u() * 0.9) / n as f64 * std::f64::consts::TAU;
            let rad = (0.05 + g.u() * 0.95) * scale;
            v.push(self.point(g, cx, cy, rad, ang));
        }
        v.dedup();
        v.push(v[0]);
        v
    }
}

fn wkb_seq(out: &mut Vec<u8>, s: &[Coord]) {
    out.extend((s.len() as u32).to_le_bytes());
    for c in s {
        out.extend(c.x.to_le_bytes());
        out.extend(c.y.to_le_bytes());
    }
}

fn wkb(geom: &Geometry) -> Vec<u8> {
    let mut o = vec![1u8];
    match geom {
        Geometry::Point(c) => {
            o.extend(1u32.to_le_bytes());
            o.extend(c.x.to_le_bytes());
            o.extend(c.y.to_le_bytes());
        }
        Geometry::LineString(s) => {
            o.extend(2u32.to_le_bytes());
            wkb_seq(&mut o, s);
        }
        Geometry::Polygon(p) => {
            o.extend(3u32.to_le_bytes());
            o.extend((p.rings.len() as u32).to_le_bytes());
            for r in &p.rings {
                wkb_seq(&mut o, r);
            }
        }
        Geometry::MultiPoint(ps) => {
            o.extend(4u32.to_le_bytes());
            o.extend((ps.len() as u32).to_le_bytes());
            for p in ps {
                o.extend(wkb(&Geometry::Point(*p)));
            }
        }
        Geometry::MultiPolygon(ps) => {
            o.extend(6u32.to_le_bytes());
            o.extend((ps.len() as u32).to_le_bytes());
            for p in ps {
                o.extend(wkb(&Geometry::Polygon(p.clone())));
            }
        }
        _ => unreachable!("not generated"),
    }
    o
}

fn hexs(b: &[u8]) -> String {
    b.iter().fold(String::with_capacity(b.len() * 2), |mut s, x| {
        let _ = write!(s, "{x:02x}");
        s
    })
}

fn cases() -> Vec<(String, GridSpec, Vec<Geometry>)> {
    let mut g = Lcg(20261031);
    let res = [0.1, 0.05, 0.3, 1.0 / 3.0, 1.0 / 7.0, 0.7, 1.0, 2.5];
    let origin = [(0.0, 0.0), (-180.0, 90.0), (0.3, 0.3), (-1.7, 2.9), (1000.1, -2000.7)];
    let mut out = Vec::new();
    for k in 0..240 {
        let (x0, y1) = g.pick(&origin);
        let grid = G { x0, y1, r: g.pick(&res), nc: 4 + (g.u() * 14.0) as u32, nr: 4 + (g.u() * 14.0) as u32 };
        let (w, h) = (grid.nc as f64 * grid.r, grid.nr as f64 * grid.r);
        let cx = x0 + (0.05 + 0.9 * g.u()) * w;
        let cy = y1 - h + (0.05 + 0.9 * g.u()) * h;
        let kind = k % 6;
        let geoms = match kind {
            // a ring, one or two of them
            0 | 1 => {
                let mut v = vec![Geometry::Polygon(Polygon::new(vec![grid.ring(&mut g, cx, cy, 1.0)]))];
                if kind == 1 {
                    let (cx2, cy2) = (x0 + g.u() * w, y1 - g.u() * h);
                    v.push(Geometry::Polygon(Polygon::new(vec![grid.ring(&mut g, cx2, cy2, 0.6)])));
                }
                v
            }
            // a ring with a hole (reversed, smaller), possibly touching it
            2 => {
                let outer = grid.ring(&mut g, cx, cy, 1.0);
                let mut inner = grid.ring(&mut g, cx, cy, 0.3);
                inner.reverse();
                vec![Geometry::Polygon(Polygon::new(vec![outer, inner]))]
            }
            // a two-part multipolygon
            3 => {
                let a = Polygon::new(vec![grid.ring(&mut g, cx, cy, 0.5)]);
                let (cx2, cy2) = (x0 + g.u() * w, y1 - g.u() * h);
                let b = Polygon::new(vec![grid.ring(&mut g, cx2, cy2, 0.5)]);
                vec![Geometry::MultiPolygon(vec![a, b])]
            }
            // an open line through the seams
            4 => {
                let mut l = grid.ring(&mut g, cx, cy, 1.0);
                l.pop();
                if l.len() < 2 {
                    l.push(Coord::new(cx, cy));
                }
                vec![Geometry::LineString(l)]
            }
            // points on corners, walls and centres
            _ => {
                let mut p = grid.ring(&mut g, cx, cy, 1.0);
                p.pop();
                vec![Geometry::MultiPoint(p)]
            }
        };
        out.push((format!("gen{k:03}"), grid.spec(), geoms));
    }
    out
}

fn tables(r: &BurnResult, s: &mut String) {
    for x in &r.runs {
        let _ = writeln!(s, "run {} {} {} {}", x.row, x.col_start, x.col_end, x.id);
    }
    for x in &r.edges {
        let _ = writeln!(s, "edge {} {} {:08x} {}", x.row, x.col, x.fraction.to_bits(), x.id);
    }
    for x in &r.lines {
        let _ = writeln!(s, "line {} {} {:08x} {}", x.row, x.col, x.length.to_bits(), x.id);
    }
    for x in &r.points {
        let _ = writeln!(s, "point {} {} {}", x.row, x.col, x.id);
    }
}

const HEADER: &str = "# canonical tables from controlledburn-rs tests/canon.rs; do not edit\n";

/// Burn one case in both modes and append its block to `s`.
fn emit(s: &mut String, name: &str, grid: &GridSpec, blobs: &[Vec<u8>]) {
    let _ = writeln!(s, "case {name}");
    let _ = writeln!(
        s,
        "grid {:016x} {:016x} {:016x} {:016x} {} {}",
        grid.xmin.to_bits(),
        grid.ymin.to_bits(),
        grid.xmax.to_bits(),
        grid.ymax.to_bits(),
        grid.ncol,
        grid.nrow
    );
    for b in blobs {
        let _ = writeln!(s, "wkb {}", hexs(b));
    }
    for (label, mode) in [("coverage", BurnMode::Coverage), ("approx", BurnMode::Approx)] {
        let r = burn_wkb(blobs.iter().map(|b| b.as_slice()), grid, BurnOptions::parity(mode)).unwrap();
        assert!(r.notes.is_empty(), "{name}: {:?}", r.notes);
        let _ = writeln!(s, "mode {label}");
        tables(&r, s);
    }
    s.push_str("end\n");
}

/// Generate the cases and render the whole file. The generator uses the
/// platform's `sin`/`cos`, so it only runs when writing the file; checks
/// read the frozen inputs back from the file.
fn render() -> String {
    let mut s = String::from(HEADER);
    for (name, grid, geoms) in cases() {
        let blobs: Vec<Vec<u8>> = geoms.iter().map(wkb).collect();
        emit(&mut s, &name, &grid, &blobs);
    }
    s
}

fn unhex(h: &str) -> Vec<u8> {
    (0..h.len()).step_by(2).map(|i| u8::from_str_radix(&h[i..i + 2], 16).unwrap()).collect()
}

/// Re-burn the inputs recorded in `committed` (grid and WKB, bit exact).
fn reburn(committed: &str) -> String {
    let mut s = String::from(HEADER);
    let mut name = String::new();
    let mut grid = None;
    let mut blobs = Vec::new();
    for line in committed.lines() {
        let mut it = line.split_whitespace();
        match it.next() {
            Some("case") => name = it.next().unwrap().to_string(),
            Some("grid") => {
                let v: Vec<&str> = it.collect();
                let f = |i: usize| f64::from_bits(u64::from_str_radix(v[i], 16).unwrap());
                grid = Some(GridSpec::new(f(0), f(1), f(2), f(3), v[4].parse().unwrap(), v[5].parse().unwrap()));
            }
            Some("wkb") => blobs.push(unhex(it.next().unwrap())),
            Some("end") => {
                emit(&mut s, &name, grid.as_ref().unwrap(), &blobs);
                grid = None;
                blobs.clear();
            }
            _ => {}
        }
    }
    s
}

#[test]
fn canonical_tables_are_current() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/canon-tables.txt");
    if std::env::var_os("CANON_WRITE").is_some() {
        std::fs::write(&path, render()).unwrap();
        return;
    }
    let committed = std::fs::read_to_string(&path).expect("fixtures/canon-tables.txt (CANON_WRITE=1 to create)");
    let committed = committed.replace("\r\n", "\n");
    let now = reburn(&committed);
    if committed != now {
        let line = committed.lines().zip(now.lines()).position(|(a, b)| a != b);
        panic!(
            "burn output differs from fixtures/canon-tables.txt (first differing line: {line:?}); \
             if the change is deliberate, regenerate with CANON_WRITE=1 and port it to the C++ core"
        );
    }
}
