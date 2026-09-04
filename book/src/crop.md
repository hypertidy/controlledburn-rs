# Cropping to a tile

[Materialize, chunk by chunk](materialize.md) turns the whole result into
pixels. Often you want the opposite: burn once, then pull out many small
windows — map tiles, a study area, the cells overlapping one Zarr chunk —
without re-burning and without ever allocating the full raster.

`crop` does that on the sparse tables directly. It filters and clips the
four tables (`runs`, `edges`, `lines`, `points`) to a target window, snaps
the window outward to whole cells, and re-bases every row/column index to 1
so the cropped result stands on its own. There is no dense buffer in sight;
the cost is proportional to the number of records, not the number of pixels.

Because a `BurnResult` carries no grid metadata of its own, you pass the
`GridSpec` the result was burned into, and `crop` hands back the snapped
sub-grid alongside the cropped tables. That sub-grid feeds straight into
`materialize`, giving the burn-once, cut-tiles workflow:

```rust
# extern crate controlledburn;
use controlledburn::{burn, crop, materialize, BurnOptions, Coord, Geometry, GridSpec, MaterializeOptions};
# use controlledburn::Polygon;

let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
let square = Geometry::Polygon(Polygon::new(vec![vec![
    Coord::new(2.5, 4.5), Coord::new(6.5, 4.5), Coord::new(6.5, 8.5), Coord::new(2.5, 8.5), Coord::new(2.5, 4.5),
]]));
let r = burn(&[square], &grid, BurnOptions::default()).unwrap();

// Cut one tile: the window x in [3, 7], y in [5, 9].
let (tile, subgrid) = crop(&r, &grid, [3.0, 5.0, 7.0, 9.0]).unwrap();
assert_eq!((subgrid.ncol, subgrid.nrow), (4, 4));

// Materialize just that tile.
let mut pixels = vec![f64::NAN; (subgrid.ncol * subgrid.nrow) as usize];
materialize(&tile, &mut pixels, subgrid.ncol, subgrid.nrow, None, &MaterializeOptions::default()).unwrap();
```

The `target` is `[xmin, ymin, xmax, ymax]`, matching the field order of
`GridSpec`. The window is snapped outward to cell boundaries, so a tile
always covers whole cells, and it is clamped to the grid — a window larger
than the grid returns the whole grid.

A window that misses the grid entirely returns `None` rather than an empty,
zero-dimension grid (which `GridSpec` would reject as invalid):

```rust
# extern crate controlledburn;
use controlledburn::{burn, crop, BurnOptions, Coord, Geometry, GridSpec};
# use controlledburn::Polygon;

let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
let square = Geometry::Polygon(Polygon::new(vec![vec![
    Coord::new(2.5, 4.5), Coord::new(6.5, 4.5), Coord::new(6.5, 8.5), Coord::new(2.5, 8.5), Coord::new(2.5, 4.5),
]]));
let r = burn(&[square], &grid, BurnOptions::default()).unwrap();

assert!(crop(&r, &grid, [100.0, 100.0, 200.0, 200.0]).is_none());
```

A crop followed by materialize is bit-for-bit the same as materializing the
full grid and slicing out the same cells — cropping only decides *which*
records survive and renumbers them, never *how* a cell is painted. Interior
`runs` that straddle the window edge are clipped to it; `edges`, `lines`, and
`points` outside the window are dropped; and `notes` are carried through
unchanged, since they are keyed by input geometry rather than by cell.

This mirrors the R package's `crop_burn()` and the Python
`BurnResult.crop()`, with one deliberate difference: those bindings order
`target` as `[xmin, xmax, ymin, ymax]`, while the Rust crate uses
`[xmin, ymin, xmax, ymax]` to line up with `GridSpec`.
