//! How closely do per-cell fractions of a polygon partition sum to 1
//! when each polygon is burned on its own?
//!
//! Builds a jittered lattice of quads that tile a square exactly (shared
//! edges have identical coordinates), burns every quad in coverage mode
//! and sums the fractions in each cell. A sum other than 1 is a gap or an
//! overlap introduced by burning shared edges twice, once in each
//! direction. See design/topology-mesh.md.
//!
//! Run with `cargo run --release --example partition_sum`.

use std::collections::HashMap;

use controlledburn::{burn, BurnOptions, Coord, Geometry, GridSpec, Polygon};

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn main() {
    let (k, dom) = (24usize, 100.0);
    let step = dom / k as f64;
    let mut rng = Lcg(7);
    // Lattice vertices, jittered in the interior, fixed on the frame.
    let mut v = vec![vec![Coord::new(0.0, 0.0); k + 1]; k + 1];
    for (j, row) in v.iter_mut().enumerate() {
        for (i, c) in row.iter_mut().enumerate() {
            let jit = |t: usize, r: f64| {
                if t == 0 || t == k {
                    0.0
                } else {
                    (r - 0.5) * 0.8 * step
                }
            };
            let (rx, ry) = (rng.next(), rng.next());
            *c = Coord::new(i as f64 * step + jit(i, rx), j as f64 * step + jit(j, ry));
        }
    }
    let mut quads = vec![];
    for j in 0..k {
        for i in 0..k {
            let ring = vec![v[j][i], v[j][i + 1], v[j + 1][i + 1], v[j + 1][i], v[j][i]];
            quads.push(Geometry::Polygon(Polygon::new(vec![ring])));
        }
    }

    for &n in &[50u32, 173, 512] {
        let grid = GridSpec::new(0.0, 0.0, dom, dom, n, n);
        let r = burn(&quads, &grid, BurnOptions::coverage()).unwrap();
        let mut sum: HashMap<(i32, i32), f64> = HashMap::new();
        for run in &r.runs {
            for col in run.col_start..run.col_end {
                *sum.entry((run.row, col)).or_default() += 1.0;
            }
        }
        for e in &r.edges {
            *sum.entry((e.row, e.col)).or_default() += e.fraction as f64;
        }
        let cells = (n as usize) * (n as usize);
        let worst = sum.values().map(|s| (s - 1.0).abs()).fold(0.0, f64::max);
        let off = sum.values().filter(|s| (*s - 1.0).abs() > 1e-6).count();
        println!(
            "grid {n} x {n}: {} of {cells} cells reached, {} edge records, max |sum - 1| = {worst:.2e}, cells off by > 1e-6: {off}",
            sum.len(),
            r.edges.len()
        );
    }
}
