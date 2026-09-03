use controlledburn::{burn, BurnOptions, Coord, Geometry, GridSpec, Polygon};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

/// Deterministic jittered "coastline" ring, same construction as the C++
/// probe used for the baseline numbers in the scoping document.
fn coastline(n: usize) -> Geometry {
    let mut state: u32 = 1;
    let mut rnd = || {
        // xorshift32, uniform in [-0.4, 0.4)
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        (state as f64 / u32::MAX as f64) * 0.8 - 0.4
    };
    let mut ring: Vec<Coord> = (0..n)
        .map(|i| {
            let a = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
            let r = 400.0 + 50.0 * (13.0 * a).sin() + rnd() * 20.0;
            Coord::new(500.0 + r * a.cos(), 500.0 + r * a.sin())
        })
        .collect();
    ring.push(ring[0]);
    Geometry::Polygon(Polygon::new(vec![ring]))
}

fn bench(c: &mut Criterion) {
    let g = coastline(200_000);
    let mut group = c.benchmark_group("coastline_200k");
    group.sample_size(10);
    for n in [1024u32, 4096] {
        let grid = GridSpec::new(0.0, 0.0, 1000.0, 1000.0, n, n);
        group.bench_with_input(BenchmarkId::new("coverage", n), &grid, |b, grid| {
            b.iter(|| burn(std::slice::from_ref(&g), grid, BurnOptions::coverage()).unwrap())
        });
        group.bench_with_input(BenchmarkId::new("approx", n), &grid, |b, grid| {
            b.iter(|| burn(std::slice::from_ref(&g), grid, BurnOptions::approx()).unwrap())
        });
    }
    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
