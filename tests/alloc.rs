//! Burning many small polygons in coverage mode reuses its buffers: past
//! the first few polygons, the only allocations are the output tables
//! growing. Before the buffers were reused this was about 90 allocations
//! per polygon here (some hundreds on real data), and about half the time
//! of a burn of many small polygons.
//!
//! A single large ring allocates per row of the sweep, not per cell it
//! crosses.
//!
//! One test in this file, as the counting allocator is global to the
//! test binary and tests run in parallel.

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};

use controlledburn::{burn, BurnOptions, Coord, Geometry, GridSpec, Polygon};

static ALLOCS: AtomicU64 = AtomicU64::new(0);

struct Counting;

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        ALLOCS.fetch_add(1, Relaxed);
        unsafe { System.alloc(l) }
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        unsafe { System.dealloc(p, l) }
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, n: usize) -> *mut u8 {
        ALLOCS.fetch_add(1, Relaxed);
        unsafe { System.realloc(p, l, n) }
    }
}

#[global_allocator]
static GLOBAL: Counting = Counting;

/// A diamond of radius `r` around `(cx, cy)`, CCW or CW, closed.
fn diamond(cx: f64, cy: f64, r: f64, ccw: bool) -> Vec<Coord> {
    let mut v = vec![Coord::new(cx, cy - r), Coord::new(cx + r, cy), Coord::new(cx, cy + r), Coord::new(cx - r, cy)];
    if !ccw {
        v.reverse();
    }
    v.push(v[0]);
    v
}

/// `n` diamonds of 1.3 to 3.7 cells across, at irrational offsets so they
/// cut cells in general position; every third has a hole.
fn polygons(n: usize) -> Vec<Geometry> {
    (0..n)
        .map(|k| {
            let (cx, cy) = (5.0 + (k % 40) as f64 * 10.0 + 0.37, 5.0 + (k / 40) as f64 * 10.0 + 0.61);
            let r = 1.3 + (k % 5) as f64 * 0.6;
            let mut rings = vec![diamond(cx, cy, r, k % 2 == 0)];
            if k % 3 == 0 {
                rings.push(diamond(cx + 0.1, cy, r * 0.4, k % 2 == 1));
            }
            Geometry::Polygon(Polygon::new(rings))
        })
        .collect()
}

fn allocs_to_burn(geoms: &[Geometry], grid: &GridSpec) -> u64 {
    let before = ALLOCS.load(Relaxed);
    let r = burn(geoms, grid, BurnOptions::coverage()).unwrap();
    let n = ALLOCS.load(Relaxed) - before;
    assert!(!r.edges.is_empty() && r.notes.is_empty());
    n
}

#[test]
fn coverage_burn_does_not_allocate_per_polygon() {
    let grid = GridSpec::new(0.0, 0.0, 400.0, 400.0, 400, 400);
    let small = polygons(400);
    let large = polygons(1600);
    let a = allocs_to_burn(&small, &grid);
    let b = allocs_to_burn(&large, &grid);
    // 1200 more polygons: only the output tables' amortised growth (a
    // handful of reallocations; 4 when written), not one or more per
    // polygon.
    let per_polygon = (b.saturating_sub(a)) as f64 / 1200.0;
    assert!(per_polygon < 0.05, "{a} allocations for 400 polygons, {b} for 1600: {per_polygon:.2} per polygon");

    // One large ring: the walk keeps its traversals in one array, so it
    // allocates per growth of that array, not per cell it crosses. What
    // is left grows with the ring's rows (the sweep's per-row tables). A
    // wobbly ring of 20,000 vertices on a fine grid crosses about 19,000
    // boundary cells in about 3,100 rows; before the traversals were flattened it
    // allocated once or more per cell.
    let fine = GridSpec::new(0.0, 0.0, 400.0, 400.0, 4000, 4000);
    let ring = |k: usize| -> Geometry {
        let mut v: Vec<Coord> = (0..k)
            .map(|i| {
                let a = i as f64 / k as f64 * std::f64::consts::TAU;
                let r = 150.0 + 7.3 * (a * 37.0).sin();
                Coord::new(200.0 + r * a.cos(), 200.0 + r * a.sin())
            })
            .collect();
        v.push(v[0]);
        Geometry::Polygon(Polygon::new(vec![v]))
    };
    let before = ALLOCS.load(Relaxed);
    let r = burn(&[ring(20_000)], &fine, BurnOptions::coverage()).unwrap();
    let c = ALLOCS.load(Relaxed) - before;
    let rows = 1 + (r.edges.iter().map(|e| e.row).max().unwrap() - r.edges.iter().map(|e| e.row).min().unwrap()) as u64;
    let cells = r.edges.len() as u64;
    assert!(cells > 15_000, "{cells} boundary cells");
    assert!(c < 5 * rows, "{c} allocations to burn one ring over {rows} rows and {cells} boundary cells");
}
