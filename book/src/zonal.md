# Zonal statistics from the tables

Zonal statistics ask, for each polygon, a summary of the raster values
under it. Done on a dense mask this is a full pass over the raster per
zone, or a rasterized zone-id layer that cannot represent shared
boundary cells honestly. Done on the tables it is a loop over `runs` and
`edges`, with the coverage fraction as the weight, and shared boundary
cells contribute their correct fraction to every zone that touches them.

## A raster to summarise

There is no raster reader in this crate, and none is needed for the
idea, so the raster here is a function of the cell: the latitude of the
cell centre on a 0.1-degree global grid. Any real raster on the same
grid is a lookup by `(row, col)` in exactly the same place.

```rust
# extern crate controlledburn;
use controlledburn::GridSpec;

let grid = GridSpec::new(-180.0, -90.0, 180.0, 90.0, 3600, 1800);

// Row 0 is the top row, so its centre is at ymax - dy / 2.
fn value(grid: &GridSpec, row: i32, _col: i32) -> f64 {
    grid.ymax - (row as f64 + 0.5) * grid.dy()
}
assert!((value(&grid, 0, 0) - 89.95).abs() < 1e-9);
assert!((value(&grid, 1799, 0) + 89.95).abs() < 1e-9);
```

## Area-weighted mean per zone

Each run contributes every cell in its range with weight 1; each edge
contributes one cell with weight `fraction`. Accumulating `sum(w * v)`
and `sum(w)` per `id` gives the mean.

```rust
# extern crate controlledburn;
use controlledburn::{burn_wkb, BurnOptions, BurnResult, GridSpec};
# let data_dir = env!("CONTROLLEDBURN_BOOK_DATA");
# let bytes = std::fs::read(format!("{data_dir}/ne_110m_countries.wkb")).unwrap();
# fn split_blobs(bytes: &[u8]) -> Vec<&[u8]> { let n = u32::from_le_bytes(bytes[0..4].try_into().unwrap()) as usize; let mut pos = 4; let mut out = Vec::with_capacity(n); for _ in 0..n { let len = u32::from_le_bytes(bytes[pos..pos + 4].try_into().unwrap()) as usize; pos += 4; out.push(&bytes[pos..pos + len]); pos += len; } out }
# let countries = split_blobs(&bytes);
# let grid = GridSpec::new(-180.0, -90.0, 180.0, 90.0, 3600, 1800);
# fn value(grid: &GridSpec, row: i32, _col: i32) -> f64 { grid.ymax - (row as f64 + 0.5) * grid.dy() }

fn weighted_mean(r: &BurnResult, n_zones: usize, grid: &GridSpec, value: impl Fn(&GridSpec, i32, i32) -> f64) -> Vec<f64> {
    let mut sum_wv = vec![0.0f64; n_zones];
    let mut sum_w = vec![0.0f64; n_zones];
    for run in &r.runs {
        let k = run.id as usize;
        for col in run.col_start..run.col_end {
            sum_wv[k] += value(grid, run.row, col);
            sum_w[k] += 1.0;
        }
    }
    for e in &r.edges {
        let k = e.id as usize;
        let w = e.fraction as f64;
        sum_wv[k] += w * value(grid, e.row, e.col);
        sum_w[k] += w;
    }
    sum_wv.iter().zip(&sum_w).map(|(a, b)| if *b > 0.0 { a / b } else { f64::NAN }).collect()
}

let r = burn_wkb(countries.iter().copied(), &grid, BurnOptions::default()).unwrap();
let mean_lat = weighted_mean(&r, countries.len(), &grid, value);

// data/ne_110m_countries.csv maps id to name; a few known ids:
let (australia, chile, iceland) = (137, 10, 144);
# // the CSV has a header line, so line `id + 1` is the row for that id
# let names = std::fs::read_to_string(format!("{data_dir}/ne_110m_countries.csv")).unwrap();
# let name_of = |id: usize| names.lines().nth(id + 1).unwrap().split(',').nth(1).unwrap().trim_matches('"').to_string();
# assert_eq!(name_of(australia), "Australia");
# assert_eq!(name_of(chile), "Chile");
# assert_eq!(name_of(iceland), "Iceland");
println!("Australia {:.2}  Chile {:.2}  Iceland {:.2}", mean_lat[australia], mean_lat[chile], mean_lat[iceland]);
# assert!((mean_lat[australia] - -25.7).abs() < 0.5, "{}", mean_lat[australia]);
# assert!((mean_lat[chile] - -39.0).abs() < 0.5, "{}", mean_lat[chile]);
# assert!((mean_lat[iceland] - 65.0).abs() < 0.5, "{}", mean_lat[iceland]);
```

Australia's area-weighted mean latitude comes out near -25.7, Chile's
near -39.0 (a long thin country weighted by its area, not its length),
Iceland's near 65.1.

## Why the fraction matters

A cell on the border between two countries is shared. With a dense
zone-id raster it belongs to one of them, and which one depends on the
rasterizer's tie rule. With the tables it belongs to both, in
proportion:

```rust
# extern crate controlledburn;
# use controlledburn::{burn_wkb, BurnOptions, GridSpec};
# let data_dir = env!("CONTROLLEDBURN_BOOK_DATA");
# let bytes = std::fs::read(format!("{data_dir}/ne_110m_countries.wkb")).unwrap();
# fn split_blobs(bytes: &[u8]) -> Vec<&[u8]> { let n = u32::from_le_bytes(bytes[0..4].try_into().unwrap()) as usize; let mut pos = 4; let mut out = Vec::with_capacity(n); for _ in 0..n { let len = u32::from_le_bytes(bytes[pos..pos + 4].try_into().unwrap()) as usize; pos += 4; out.push(&bytes[pos..pos + len]); pos += len; } out }
# let countries = split_blobs(&bytes);
# let grid = GridSpec::new(-180.0, -90.0, 180.0, 90.0, 3600, 1800);
# let r = burn_wkb(countries.iter().copied(), &grid, BurnOptions::default()).unwrap();
use std::collections::HashMap;

// Total coverage per cell across all countries, from edges only.
let mut per_cell: HashMap<(i32, i32), (f64, u32)> = HashMap::new();
for e in &r.edges {
    let c = per_cell.entry((e.row, e.col)).or_insert((0.0, 0));
    c.0 += e.fraction as f64;
    c.1 += 1;
}
let shared = per_cell.values().filter(|(_, n)| *n > 1).count();
let over_full = per_cell.values().filter(|(f, _)| *f > 1.0 + 1e-4).count();
println!("{} boundary cells are shared by two or more countries", shared);
# assert!(shared > 10_000);
// Neighbours that share a border tile the cell: their fractions sum to at most 1
// (tiny excesses come from Natural Earth's own boundary slivers).
assert!(over_full < shared / 100);
```

About 24,000 cells on this grid are shared, and in all but four of them
the neighbouring countries' fractions sum to at most one (exactly one
where the border is the only thing in the cell): the cell is divided,
not assigned. That is the difference between a coverage
representation and a mask, and it is why the same tables serve both
rasterization and exact zonal statistics.
