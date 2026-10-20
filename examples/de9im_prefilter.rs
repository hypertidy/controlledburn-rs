//! Cell-level prefilter for polygon-polygon predicates on a shared grid.
//!
//! Burns a set of random polygons in coverage mode and uses only the runs
//! (fully covered cells) and edges (partially covered cells) to answer
//! `intersects`, "interiors meet" and `within` for every bbox-overlapping
//! pair as true, false or undecided. Decided answers are checked against
//! `geo::Relate`. See design/de9im-cells.md.
//!
//! Run with `cargo run --release --example de9im_prefilter`.

use std::collections::{HashMap, HashSet};

use controlledburn::{burn, BurnOptions, Coord, Geometry, GridSpec, Polygon};
use geo::Relate;

/// Small deterministic generator, so the example needs no extra crates.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// A random star-shaped ring inside [0, dom] x [0, dom].
fn star(rng: &mut Lcg, dom: f64) -> Vec<(f64, f64)> {
    let cx = 5.0 + rng.next() * (dom - 10.0);
    let cy = 5.0 + rng.next() * (dom - 10.0);
    let rad = 1.0 + rng.next() * 12.0;
    let n = 8 + (rng.next() * 40.0) as usize;
    let mut ring: Vec<(f64, f64)> = (0..n)
        .map(|i| {
            let t = i as f64 / n as f64 * std::f64::consts::TAU;
            let r = rad * (0.4 + 0.6 * rng.next());
            ((cx + r * t.cos()).clamp(0.0, dom), (cy + r * t.sin()).clamp(0.0, dom))
        })
        .collect();
    ring.push(ring[0]);
    ring
}

/// Cells of one polygon: fully covered (runs) and partially covered (edges).
#[derive(Default)]
struct Cells {
    full: HashSet<(i32, i32)>,
    part: HashSet<(i32, i32)>,
}

impl Cells {
    fn any(&self, k: &(i32, i32)) -> bool {
        self.full.contains(k) || self.part.contains(k)
    }

    fn all(&self) -> impl Iterator<Item = &(i32, i32)> {
        self.full.iter().chain(self.part.iter())
    }
}

/// A full cell of one polygon that the other has any area in.
fn full_meets_any(a: &Cells, b: &Cells) -> bool {
    a.full.iter().any(|k| b.any(k)) || b.full.iter().any(|k| a.any(k))
}

/// `Some(answer)` when the cells decide it, `None` when exact refinement
/// is needed. Without a touches table, "false" needs a one-cell dilation:
/// polygons meeting along a cell wall share no cell.
fn intersects(a: &Cells, b: &Cells) -> Option<bool> {
    if full_meets_any(a, b) {
        return Some(true);
    }
    let near = a.all().any(|&(r, c)| (-1..=1).any(|dr| (-1..=1).any(|dc| b.any(&(r + dr, c + dc)))));
    if near {
        None
    } else {
        Some(false)
    }
}

/// Interiors meet (intersects and not touches). Positive-area overlap
/// must show up as area of both in some cell, so no dilation is needed.
fn interiors_meet(a: &Cells, b: &Cells) -> Option<bool> {
    if full_meets_any(a, b) {
        return Some(true);
    }
    if a.part.iter().any(|k| b.part.contains(k)) {
        None
    } else {
        Some(false)
    }
}

/// A within B.
fn within(a: &Cells, b: &Cells) -> Option<bool> {
    if a.full.iter().any(|k| !b.full.contains(k)) || a.part.iter().any(|k| !b.any(k)) {
        return Some(false);
    }
    if a.all().all(|k| b.full.contains(k)) {
        return Some(true);
    }
    None
}

fn main() {
    let dom = 100.0;
    let mut rng = Lcg(42);
    let rings: Vec<Vec<(f64, f64)>> = (0..300).map(|_| star(&mut rng, dom)).collect();
    let exact_polys: Vec<geo::Polygon> = rings.iter().map(|r| geo::Polygon::new(r.clone().into(), vec![])).collect();
    let geoms: Vec<Geometry> = rings
        .iter()
        .map(|r| Geometry::Polygon(Polygon::new(vec![r.iter().map(|&(x, y)| Coord::new(x, y)).collect()])))
        .collect();

    // The baseline a spatial index already gives: bbox-overlapping pairs.
    let bbox: Vec<(f64, f64, f64, f64)> = rings
        .iter()
        .map(|r| {
            r.iter().fold((f64::MAX, f64::MAX, f64::MIN, f64::MIN), |b, &(x, y)| {
                (b.0.min(x), b.1.min(y), b.2.max(x), b.3.max(y))
            })
        })
        .collect();
    let mut pairs = vec![];
    for (i, a) in bbox.iter().enumerate() {
        for (j, b) in bbox.iter().enumerate() {
            if i != j && a.0 <= b.2 && b.0 <= a.2 && a.1 <= b.3 && b.1 <= a.3 {
                pairs.push((i, j));
            }
        }
    }

    let t = std::time::Instant::now();
    let truth: Vec<[bool; 3]> = pairs
        .iter()
        .map(|&(i, j)| {
            let m = exact_polys[i].relate(&exact_polys[j]);
            [m.is_intersects(), m.is_intersects() && !m.is_touches(), m.is_within()]
        })
        .collect();
    println!(
        "{} polygons, {} bbox-overlapping ordered pairs, geo Relate on all pairs: {:?}",
        rings.len(),
        pairs.len(),
        t.elapsed()
    );

    let names = ["intersects", "interiors meet", "within"];
    for &n in &[32u32, 128, 512] {
        let grid = GridSpec::new(0.0, 0.0, dom, dom, n, n);
        let t = std::time::Instant::now();
        let r = burn(&geoms, &grid, BurnOptions::coverage()).unwrap();
        let t_burn = t.elapsed();

        let mut cells: HashMap<i32, Cells> = HashMap::new();
        for run in &r.runs {
            let c = cells.entry(run.id).or_default();
            for col in run.col_start..run.col_end {
                c.full.insert((run.row, col));
            }
        }
        for e in &r.edges {
            cells.entry(e.id).or_default().part.insert((e.row, e.col));
        }

        // (decided, wrong) per predicate
        let mut stats = [(0usize, 0usize); 3];
        for (k, &(i, j)) in pairs.iter().enumerate() {
            let (a, b) = (&cells[&(i as i32)], &cells[&(j as i32)]);
            let got = [intersects(a, b), interiors_meet(a, b), within(a, b)];
            for p in 0..3 {
                if let Some(g) = got[p] {
                    stats[p].0 += 1;
                    if g != truth[k][p] {
                        stats[p].1 += 1;
                    }
                }
            }
        }
        println!("grid {n} x {n}: burn {t_burn:?}");
        for (name, (decided, wrong)) in names.iter().zip(stats) {
            println!("  {:16} decided {:5.1}%  wrong {}", name, 100.0 * decided as f64 / pairs.len() as f64, wrong);
            assert_eq!(wrong, 0, "cell rule for {name} disagrees with geo::Relate");
        }
    }
}
