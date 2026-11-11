# Tile-parallel burning, and routing geometries by footprint

Status: design, not implemented. Scope: how controlledburn-rs should
produce tiled output in parallel (per chunk of a data cube, per worker,
per machine), drawing on the measurements in
[hypertidy/cb-compare](https://github.com/hypertidy/cb-compare)
([writeup](https://hypertidy.github.io/cb-compare/)). Audience:
developers of controlledburn, burnscars, and front ends that write tiled
or chunked output (Zarr, COG, xarray).

## Summary

There are two ways to get tiled output from a set of geometries:

- **Burn once, then slice.** Burn every geometry over the whole grid,
  then route the records to tiles by integer arithmetic
  (`slice_tiles`, `slice_tiles_par`).
- **Burn per tile.** For each tile, burn the geometries whose bounding
  box meets it onto that tile's window of the grid.

Burning per tile with simple batching is the natural design for parallel
and lazy work, and it needs no cutting of the vectors. Its one bad case
is a large geometry over many small tiles, whose vertices are re-read
once per tile. Burning once has no such case, but holds the whole output
and doesn't parallelize a single geometry.

The proposal is to use both, choosing per geometry by its footprint:
**geometries that touch few tiles are batched and burned per tile;
geometries that touch many are burned once and sliced.** The bounding-box
index the batching needs already gives each geometry's tile count. This
keeps the simple design for the common case and removes its bad case,
without the prepared-geometry machinery that would otherwise be needed.

For polygons this works in controlledburn-rs today, with no new code:
burn each tile's batch onto the tile's own `GridSpec` (the one
`slice_tiles` reports). Checked on real data, the tiles match the
whole-grid burn record for record, bit for bit. Points and lines need
one small addition, a burn onto a *window of a grid* that keeps the
full grid's lattice (below, "Exact tile windows").

## Burning per tile

Split the output grid into tiles (the chunks of the array being
written). For each tile:

1. take the geometries whose bounding box meets the tile, from a spatial
   index built once (a list per tile, or an R-tree);
2. burn them onto the tile's window of the grid;
3. write the tile.

The vectors are never cut. The engine already clips to the grid it is
given, so a geometry that extends past the tile costs nothing outside
it, and a geometry that spans several tiles is simply burned once in
each. Each tile's output is complete and independent of every other
tile's, so tiles can be burned:

- in parallel on one machine, one tile per task;
- across machines, with no shared state but the input;
- lazily, when a chunk of a data cube is first read, with no global
  step at all.

A worker's memory is bounded by its tile's output and its geometries,
not by the whole grid. This is the design burnscars calls
*raster-parallel* (brix chunks), and it is what cb-compare measures as
"re-burn per fixed tile".

### What it rests on

**Exact tile windows.** A tile must be a window of the full grid's
lattice, not a new grid built from the tile's floating-point extent.
Every cell wall and cell centre has to come from the full grid's origin,
step and global index (`xmin + col * dx`), so that:

- a coverage fraction in a tile is the same `f32` as in the whole-grid
  burn;
- the cell-centre rule (approx mode) decides each centre the same way;
- a line is cut at the same walls, so its per-cell lengths match;
- a point or a vertex on a seam between tiles lands in exactly one tile.
  The window's own right and bottom edges are open; only the full
  grid's are closed.

How close controlledburn-rs gets today, burning each tile's batch onto
the tile's own `GridSpec` and comparing with the whole-grid burn sliced
into the same tiles:

- **Polygons: identical.** Every tile and every record, in both modes,
  on cb-compare's `rects` (1,997 rectangles on a non-dyadic 1000 × 700
  grid, a quarter with corners on lattice lines), and the rock outcrops
  at 1024 m and 256 m (25,954 polygons); about 1.8 M records in all.
- **Points on a seam: counted twice.** A tile's own grid closes its own
  right and bottom edges, so a point exactly on a seam lands in both
  tiles (37 extra points on the 1024 m rock grid, with points on every
  seventh lattice node). A caller can fix this today: drop a tile's
  points on its right and bottom edges, unless that edge is the grid's.
- **Points and lines exactly on an interior wall: can move to the
  neighbouring cell** when the cell size is not exact in binary floating
  point. The tile recomputes its own step from its extent, which differs
  from the full grid's in the last bits, so a line lying on a wall falls
  into the row on the other side of it (64-300 line records on `rects`,
  even in the top-left tile). A caller can't fix this; it needs the
  window burn below.

burnscars handles both with `Grid::window()`, which keeps the full
grid's lattice and edge rules in every window.

**Cost follows vertices × tiles touched.** Each tile re-reads every
vertex of every geometry that meets it. For many small geometries that
is close to free: most touch one tile, a few touch two or four. For a
large geometry over small tiles it is the whole cost.

### The window burn (the controlledburn-rs piece)

A sketch, not a settled API:

```rust
/// Burn `geoms` onto the window `rows` x `cols` of `grid`: the same
/// records, with the same values, as the whole-grid burn has in those
/// cells. Indices are the full grid's, so windows concatenate.
pub fn burn_window(
    geoms: &[Geometry],
    grid: &GridSpec,
    rows: Range<u32>,
    cols: Range<u32>,
    opts: BurnOptions,
) -> Result<BurnResult, BurnError>;
```

The acceptance test is a property test: for random geometries, grids and
tilings, the tiles' records put together equal the whole-grid burn, bit
for bit, in both modes, for polygons, lines and points (with points on
seams, lines along walls, and vertices on lattice lines among the
cases). burnscars tests the same property for its windows, and
controlledburn-rs's `slice_tiles` is already tested against the
whole-grid result the same way (`tests/tiles.rs`).

Internally the engine already burns onto a sub-grid snapped to the full
grid (`Grid::shrink_to_fit`, from exactextract), so the work is in
exposing a window whose edges follow the full grid's rules, not in the
walker or the sweep. For polygons it changes nothing (they already
match); it is for points and lines on seams and walls.

## What the measurements say

From cb-compare's run on one Setonix node (2 × AMD EPYC 7763, 128 cores;
Slurm job 52790194, controlledburn-rs `e1696fc`, burnscars `1d24ef7`).
"Per tile" is burnscars' re-burn per fixed tile; "once" is
controlledburn-rs's whole-grid burn plus `slice_tiles_par`. Times are
medians of 5.

**Many small geometries, a whole node: per tile wins.**

| input | mode | tiles | per tile, 128 threads | once, 128 threads |
|---|---|---|---|---|
| rock outcrops, 32 m grid (25,954 polygons) | coverage | 4,138 (512-cell) | 0.036 s | 0.13 s |
| small-1e5 (100,000 small polygons) | coverage | 1,024 (256-cell) | 0.065 s | 0.16 s |
| rock outcrops, 32 m grid | approx | 640 (2048-cell) | 0.030 s | 0.048 s |

Per tile also used far less memory: 44 MB at 1 thread on the rock
outcrops at 32 m, against 0.35 GB to hold the whole-grid output.

**One large geometry: it depends on the tile size.** coast-1e6 is a
single ring of 10⁶ vertices, in coverage mode:

| tile | tiles touched | per tile, 1 thread | per tile, 128 threads | once, 128 threads | per-tile peak memory, 128 threads |
|---|---|---|---|---|---|
| 4096 | 16 | 33.8 s | 4.1 s | 12.8 s | 45 GB |
| 1024 | 178 | 117 s | 8.9 s | 13.0 s | 73 GB |
| 256 | 2,317 | 463 s | 32.0 s | 12.9 s | 35 GB |

Per tile is the only strategy that uses the cores on a single geometry,
so with large tiles it wins by 3×. With small tiles the cost of reading
10⁶ vertices 2,317 times dominates, and burning once wins at every tile
size on one thread (2.6× to 37×). Burning once uses 13.7 GB here.

**Few tiles cap the parallelism.** With 4096-cell tiles, small-1e5 has
only 4 tiles, so per tile can use only 4 workers: 0.47 s at 128
threads, against 0.23 s to burn once.

## Routing by footprint

Both strategies produce records keyed by `(id, row, col)` in the full
grid's indices. So a tile's output can come from both at once, with no
merge, because no geometry goes through both:

1. Build the spatial index: for each geometry, the tiles its bounding
   box meets.
2. A geometry that meets at most `K` tiles goes to those tiles' batches.
3. The rest are burned once over the whole grid (feature-parallel, with
   `id_offset` for workers' slices), and the result is sliced
   (`slice_tiles_par`).
4. Each tile's output is its own batch's burn plus its slice of the
   large geometries' burn: concatenated, then sorted by id if the writer
   wants id order.

`K` is a cost trade-off: per tile costs about vertices × tiles touched,
and burning once costs about vertices plus the slice. It should be
measured, not guessed. The data so far suggest it is small (a geometry
touching 16 tiles was still 3× faster per tile on 128 threads, while
178 or more was slower). It could also be stated in vertices × tiles
rather than tiles alone, since a small polygon that touches four tiles
is cheap either way.

For lazy, chunk-on-demand use, step 3 is the one global step: the large
geometries' burn is done once and cached, and each chunk then takes its
slice. If there are no large geometries, there is no global step.

The bounding box overestimates the tiles a geometry touches, such as a
long diagonal line or an L-shaped polygon. That costs only some empty
burns, never a wrong answer.

## Not covered here

- **Prepared geometry**, as in burnscars' plans R3 and R4 (share the
  prepared edges of a large geometry across tiles, and pre-clip each
  ring to the tile). This is the other way to remove per tile's bad
  case. Routing by footprint avoids needing it. The two can coexist.
- **Grids that aren't regular lattices.** Tiles of a rotated,
  rectilinear or curvilinear grid are tiles of its index space
  (`design/grid-frontends.md`), so everything here applies after the
  front end's transform. A front end that folds columns back (periodic
  or tripolar grids) re-merges records after the fold, by the rule
  `burn` uses for multi-part geometries (`merge_cells`, in the pending
  API changes for 0.3.0).
- **Writing the tiles** (Zarr, COG), and how chunk boundaries line up
  with the tiles. A tile here can be any rectangle of whole cells.
