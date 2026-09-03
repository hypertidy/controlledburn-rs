# Materialize, chunk by chunk

*Coming in a later revision.* The tables can be turned into pixels one
window at a time, which is how a chunked format such as Zarr is
written without ever holding the full raster: the 0.01-degree world of
the [sparse versus dense](sparse-vs-dense.md) chapter is 5 GB dense
and never needs to exist in memory.

```rust
# extern crate controlledburn;
use controlledburn::{burn, materialize, BurnOptions, Coord, Geometry, GridSpec, MaterializeOptions};

let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
let square = Geometry::Polygon(Polygon::new(vec![vec![
    Coord::new(2.5, 4.5), Coord::new(6.5, 4.5), Coord::new(6.5, 8.5), Coord::new(2.5, 8.5), Coord::new(2.5, 4.5),
]]));
# use controlledburn::Polygon;
let r = burn(&[square], &grid, BurnOptions::default()).unwrap();
let mut pixels = vec![f64::NAN; 100];
materialize(&r, &mut pixels, 10, 10, None, &MaterializeOptions::default()).unwrap();
assert_eq!(pixels.iter().filter(|v| !v.is_nan()).count(), 21);
```
