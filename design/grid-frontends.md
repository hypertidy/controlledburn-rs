# Coverage on grids that are not a regular lattice

Status: discussion. Audience: developers evaluating controlledburn as the
polygon-to-cell engine for array libraries (xarray, rasterix, Zarr-backed
data cubes), where grids are often rotated, rectilinear or curvilinear
rather than a regular lattice.

## Summary

controlledburn burns polygons onto a regular lattice and returns sparse
runs of whole cells plus exact coverage fractions for boundary cells. It
does not need to know about any other kind of grid. A small front end,
outside the engine, can move each polygon into the grid's **index space**
(fractional column, row) and burn it on a lattice of unit cells; the
results are then already in the original grid's `(row, col)`.

How exact that is depends only on whether the grid is affine within each
cell:

| grid | front end | result |
|---|---|---|
| regular | none | exact (the engine's own case) |
| rotated or sheared affine | inverse geotransform, vertices only | exact |
| rectilinear (1D coordinates, uneven spacing) | cut edges at cell walls, map each axis | exact |
| curvilinear (2D coordinates) | cut edges into short pieces, inverse lookup | close, not exact |
| any of the above, polygon in another CRS | reproject, then as above, with edges cut | as above |

"Exact" here means the fractions match a world-space reference (each cell
outline clipped by the polygon in the polygon's CRS) to the engine's
single precision. For curvilinear grids no front end makes index-space
fractions exact, because a fraction of a unit cell is not a fraction of
the cell's area when the map is not affine within the cell. The error is
around one percent of a cell for a typical polar sector and is measured
below. An exact answer there needs the world-space clip.

No change to the engine is needed for any of this. A native rectilinear
axis in the engine (the grid-geometry seam) would remain worthwhile for
speed and precision, not for correctness.

## Index space as the common currency

The engine's contract is a regular grid: an extent and a column and row
count, row 0 at the top. Index space is such a grid. With `w` columns and
`h` rows:

```rust
use controlledburn::{burn, BurnOptions, Coord, Geometry, GridSpec, Polygon};

// unit cells, x = col, y = -row (so row 0 is at the top, as the engine expects)
let grid = GridSpec::new(0.0, -(h as f64), w as f64, 0.0, w as u32, h as u32);

// each ring: positions -> fractional (col, row) by the front end, then y negated
let rings: Vec<Vec<Coord>> = polygon
    .iter()
    .map(|r| to_index_ring(r).iter().map(|p| Coord::new(p[0], -p[1])).collect())
    .collect();

let r = burn(&[Geometry::Polygon(Polygon::new(rings))], &grid, BurnOptions::coverage())?;
// r.runs and r.edges carry (row, col) of the original grid, unchanged
```

Everything specific to the grid lives in `to_index_ring`. It has two jobs:

1. map positions to fractional `(col, row)`;
2. cut edges where needed, so that each piece is a straight segment in
   index space (a straight edge in the polygon's CRS is generally a curve
   in index space).

![The same polygons in a polar sector's CRS and in its index space, with vertices moved alone and with edges cut](img/two-spaces.png)

On the left, a curvilinear grid (a polar sector) in its own CRS, cells
shaded by covered fraction. In the middle, the polygons with only their
vertices moved to index space: the long sliver is moved end to end in a
single step and cuts straight across cells it never touches. On the
right, the edges are cut into short pieces first and follow the grid.
Moving vertices alone (what GDAL's rasterizer does when the target is
not a regular lattice in the polygon's CRS) is wrong by whole cells along
long edges on any non-affine grid.

## Grid families

### Rotated or sheared affine

A geotransform `[x0, dx, rx, y0, ry, dy]` with nonzero `rx` or `ry`, as
rasterix's `RasterIndex` can hold. Its inverse maps positions to index
space. An affine map keeps straight lines straight and scales every area
by the same factor, so moving the vertices alone is exact, and the
fractions are fractions of each cell's area. Nothing else is needed.

### Rectilinear

1D coordinates along x and y with uneven spacing: a stretched model grid,
or any `DataArray` whose dimension coordinates are not evenly spaced. The
front end needs cell **walls**. With CF `bounds` they are given. With
centres only, they must be chosen; midpoints between neighbouring centres
(and half a cell beyond each end) are the usual choice. This is a
convention, not a recovery: the true walls cannot be derived from centres
unless the spacing is even.

With walls known, each axis maps linearly within each cell, so the map is
affine within each cell. Cutting every edge where it crosses a wall makes
every piece straight in index space, and the burn is exact, area fractions
included. Each axis maps on its own (a binary search into the walls, then
a linear step), so the front end is a few dozen lines.

### Curvilinear

2D coordinate arrays: swaths, ocean and atmosphere model grids, rotated
pole grids read as 2D longitude and latitude. Two things change.

- **The inverse is a search.** Mapping a position to `(col, row)` means
  finding the cell that holds it (a backmap, a k-d tree, or a walk) and
  inverting the bilinear interpolation inside it. That belongs to the
  front end, or to the library that already holds the coordinates.
- **The measure is not area.** Within a cell the map is bilinear (or, for
  an analytic grid, whatever the projection is), not affine. Cutting edges
  into short pieces makes the polygon's shape right in index space, but the
  engine then reports the fraction of the unit square covered, which is
  not the fraction of the cell's area.

![Index-space fraction (edges cut) minus world-space fraction, per cell](img/difference.png)

The difference is zero for rotated affine and for rectilinear grids with
walls given. Given only centres, a rectilinear grid read through bilinear
interpolation of its centres has cells that straddle two interpolation
patches in each direction, and it behaves like a curvilinear one. On the
polar sector the difference is up to 1.5 percent of a cell, largest on
the inner rows where cells are most trapezoidal. Totals stay close
(about 1e-4 of the polygon's area), because the errors along a boundary
partly cancel.

![Largest per-cell difference as edges are cut finer, and the cost of each method](img/densify.png)

On the left, cutting edges finer converges for rectilinear grids with
walls given and stops improving on the others: what is left is the
measure, not the shape. On the right, cutting edges costs nothing
measurable, and the world-space clip costs a few times the burn.

An exact curvilinear answer clips each cell outline (traced from the
coordinates: corners, plus side midpoints when the coordinates are cell
centres) by the polygon in the polygon's CRS. The index-space burn is
still useful there, to find the handful of boundary cells worth clipping.

## Evidence

`examples/grid_frontends.rs` builds a rotated affine grid, a rectilinear
grid and a polar sector (with an analytic inverse, standing in for a
coordinate search), burns a star with a hole and a long sliver on each
through the front ends above, and compares every cell with a world-space
clip of the true cell outline. It depends only on this crate.

```text
cargo run --release --example grid_frontends

grid             polygon  area error, ref area error, burn  largest cell diff      vertices only
rotated affine         0          2.2e-16         -8.9e-10             0.0000             0.0000
rotated affine         1          9.8e-15          -1.8e-9             0.0000             0.0000
rectilinear            0         -4.4e-16          -4.2e-9             0.0000             0.1220
rectilinear            1         -2.3e-15           1.5e-9             0.0000             1.0000
polar sector           0          3.9e-14           1.8e-4             0.0075             0.3065
polar sector           1         -3.4e-14          -1.5e-4             0.0147             0.9954
```

- "area error, ref": the world-space clip's covered area against the
  polygon's area. It is exact, which validates the reference.
- "area error, burn": the burn's fractions weighted by each cell's world
  area, against the polygon's area.
- "largest cell diff": the largest per-cell difference in fraction
  between the burn (edges cut) and the reference.
- "vertices only": the same with vertices moved alone.

The figures above come from an independent implementation of the same
two methods (including grids given by centres and by corners), whose
world-space clip agrees with OGR/GEOS intersection areas to 1e-12 of a
cell.

## For array libraries

What each common xarray representation needs from a front end:

| representation | positions to index | edges cut | fractions |
|---|---|---|---|
| affine transform (e.g. rasterix `RasterIndex`), any rotation | inverse affine | no | exact |
| 1D dimension coordinates with CF bounds | per-axis search into walls | at walls | exact |
| 1D dimension coordinates, centres only | as above, walls at midpoints | at walls | exact for the chosen walls |
| 2D coordinates (curvilinear) | cell search plus inverse bilinear | short pieces | index-space; world clip for exact |

Because results are in the original grid's `(row, col)`, they index
straight into the array (`isel` positions) and into its chunks: a run
`col_start..col_end` on row `r` falls into chunks by integer division,
with no reprojection or resampling of the data. Sparse output means a
zone's weights are the size of its boundary, not of the grid.

A polygon in a different CRS composes with all of this: reproject to the
grid's CRS and map to index space as one function, and cut edges, since a
reprojected straight edge is itself a curve.

## For the engine

1. **No change is required.** The front ends are a few dozen lines each
   and sit naturally in the binding layer (Python, R) or in the array
   library, which already holds the coordinates.
2. **Helpers could ship with the crate**, without touching the burn: for
   example `burn_affine(geoms, geotransform, shape, opts)` and
   `burn_rectilinear(geoms, x_walls, y_walls, opts)`, each a front end plus
   a call to `burn`. They would give every binding the same exact
   behaviour and keep the index-space trick out of user code.
3. **A native rectilinear axis** (cell walls as a list, in place of
   `xmin, dx`) would let the engine burn rectilinear grids in their own
   CRS: no transform pass, no cutting at walls, and no loss of precision
   from moving coordinates into index space. This is an optimisation and
   a precision gain, not a correctness fix.
4. **Curvilinear stays out of the engine.** The inverse needs the
   coordinate arrays and a search structure, and exact fractions need a
   world-space clip. The engine's role is the fast first pass: which
   cells, and approximately how much.

## Open questions

- Fractions are single precision. For exact zonal sums on very large
  cells this may matter more than the curvilinear measure does.
- Which helpers, if any, belong in this crate rather than in the
  bindings.
- For curvilinear grids, whether index-space fractions are acceptable as
  documented weights (they are consistent, and close in total), or
  whether a world-space clip should be offered as an option.
- Periodic and wrapped grids (global longitude, tripolar folds) need the
  front end to split polygons at the seam; this is the same problem the
  engine has on a regular global grid.
