# Lines and points

*Coming in a later revision.* Lines produce one record per cell with the
length of line inside it; points produce one record per point. The
schema is the same shape as for polygons and the burn mode does not
affect either.

```rust
# extern crate controlledburn;
use controlledburn::{burn, BurnOptions, Coord, Geometry, GridSpec};

let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
let line = Geometry::LineString(vec![Coord::new(0.5, 0.5), Coord::new(9.5, 7.5)]);
let r = burn(&[line], &grid, BurnOptions::default()).unwrap();
let total: f64 = r.lines.iter().map(|l| l.length as f64).sum();
assert!((total - (81.0f64 + 49.0).sqrt()).abs() < 1e-4);
```
