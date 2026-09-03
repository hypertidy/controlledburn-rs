# Coverage versus approx

*Coming in a later revision.* The crate offers two answers to the
boundary question: `BurnMode::Coverage` (exact area fractions, the
default) and `BurnMode::Approx` (a cell is inside iff its centre is,
the rule GDAL and fasterize use). This chapter will put the two side by
side on a coastline and show where and why they disagree.

```rust
# extern crate controlledburn;
use controlledburn::{burn, BurnOptions, Coord, Geometry, GridSpec, Polygon};

let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
let square = Geometry::Polygon(Polygon::new(vec![vec![
    Coord::new(2.5, 4.5), Coord::new(6.5, 4.5), Coord::new(6.5, 8.5), Coord::new(2.5, 8.5), Coord::new(2.5, 4.5),
]]));
let exact = burn(std::slice::from_ref(&square), &grid, BurnOptions::coverage()).unwrap();
let approx = burn(std::slice::from_ref(&square), &grid, BurnOptions::approx()).unwrap();
assert_eq!(exact.covered_cells(), 16.0);   // 9 full + 16 partial cells
assert_eq!(approx.run_cells(), 16);        // 16 full cells, no edges
assert!(approx.edges.is_empty());
```
