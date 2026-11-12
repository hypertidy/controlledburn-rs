# Changelog

## Unreleased

- **Breaking: `BurnOptions`, `BurnError` and `WkbError` are
  `#[non_exhaustive]`**, so later versions can add options and error
  variants without breaking callers. Build options with `default()`,
  `coverage()`, `approx()`, `parity()` and the `with_*` methods (the
  fields can still be read and set); a `match` on either error needs a
  `_` arm.
- **`burn_window`: a burn onto a window of the grid**, rows by columns,
  with the records `burn` gives in those cells, the same values in the
  same order and the full grid's indices, so each tile of a tiling can
  burn just the geometries that meet it (`design/tile-parallel.md`). It
  keeps the full grid's lattice and edge rules: a point on a seam lands
  in one window, and a line on a wall stays in the whole burn's cell,
  which burning onto each tile's own `GridSpec` does not give. Exact by
  construction: a polygon is burned on the sub-grid the whole burn uses,
  sweeping only the window's rows, and a line is walked along its parts
  near the window, each started at a vertex clear of every wall. Two
  property tests check every window of random tilings against
  `slice_tiles` of the whole burn, record for record (both modes,
  coalescing and part merging on and off, lattice-aligned rectangles,
  lines along walls and points on corners, grid origins up to 1e7). A
  window costs about the vertices of the geometries that meet it; a
  polygon much larger than the window is walked whole for every window
  it meets (coast-1e5 over 177 tiles: 72 s, against 2 s on per-tile
  grids), so burn those once and slice.
- **`merge_cells`: one record per cell and id for a whole result**, by
  the rule `burn` uses for multi-part geometries. For records that come
  together after a burn, such as a curvilinear or periodic front end
  folding columns back (`design/grid-frontends.md`). A property test
  checks that burning with `merge_parts` off and then merging gives the
  default burn bit for bit.
- **`slice_tiles_par`: `slice_tiles` on several threads**, with the same
  tiles and records in the same order (std threads only; `threads = 0` for
  all available, which is 1 on targets without threads). Records are
  partitioned by band of tile rows in burn order, one `u32` index per
  record, and the bands are sliced in parallel by the serial code. On the
  dev VM, slicing the rock outcrops at 32 m (8.4 M records) went from
  0.12-0.20 s to 0.035-0.04 s on 16-32 threads, and small-1e5 (13.6 M
  records) from 0.36 s to 0.05 s at 256-cell tiles. A property test checks
  it against `slice_tiles` for random inputs, tile shapes and thread
  counts. `Tile` now derives `PartialEq`. It uses at most one thread per
  131,072 records (serial below 262,144): on Setonix, starting 128
  threads cost about 15 ms, so small outputs sliced up to 30x slower than
  serially. On the dev VM at 32 threads the cap makes those 2-10x faster
  (rock outcrops at 1024 m: 5.9 ms to 0.56 ms) and outputs of 0.2-2 M
  records 1.3-1.6x faster, with large outputs unchanged.
- **An external oracle, and more properties** (from burnscars, Apache-2.0;
  `fixtures/oracle/`). `tests/oracle.rs` checks burnscars' 50 cases against
  exactextract 0.3.0 and GDAL's centre rule, and requires each case's
  differing cells to be exactly the expected set: the 5 coverage cells
  where exactextract 0.3.0 gives 0 for a full interior cell, and the
  centre-on-edge ties burnscars lists, less 9 where cb agrees with GDAL.
  New properties: a ring and its complement sum to 1 in every cell under
  any orientations; a polygon with a hole conserves area and matches the
  per-cell clipping oracle.
- **One record per cell for multi-part geometries** (`BurnOptions::merge_parts`,
  on by default). The parts of a multipolygon were burned one by one, so
  a cell two parts cross had two edges: two halves either side of a shared
  border came out as 0.5 + 0.5 rather than a full cell, and an overlap
  stacked an edge on a run. The parts' edges are now summed per cell (in
  part order), a sum reaching a full cell becomes a run, edges inside
  another part's run are dropped, and the runs are coalesced. Line parts
  are summed per cell the same way. Overlapping parts are capped at a full
  cell, not unioned. On the rock outcrops at 32 m this merges 18 records,
  and the 14 cells that differed from burnscars by float rounding are now
  bit-equal. `parity()` turns it off: the C++ core keeps per-part records.
- **Notes for skipped rings and lines.** A ring of fewer than 4
  coordinates once closed (its part, if it is an exterior; itself, if a
  hole) and a one-coordinate line were skipped silently; each now gets a
  note, and the rest of the geometry is still burned. An all-NaN point
  (WKB `POINT EMPTY`) is empty: it was noted as non-finite.
- **`BurnOptions::id_offset`.** Ids are `id_offset + k`, notes included, so
  burns of consecutive slices of a collection concatenate without
  renumbering. A negative offset, or ids past `i32::MAX`, is the new
  `BurnError::InvalidOptions`. `BurnOptions` gains two public fields, so
  struct literals need `..Default::default()`.
- **Coverage mode reuses its buffers from one polygon to the next**
  rather than allocating per ring, row and boundary cell: about 2x faster
  on many small polygons (small-1e5 4.11 s to 1.96 s, rock outcrops at
  32 m 2.27 s to 1.20 s), level with burnscars, output bit-identical.
  `tests/alloc.rs` pins it.
- **`slice_tiles`: a burn cut into a regular tile grid in one pass.**
  Each record goes to tile `(row / tile_nrow, col / tile_ncol)` with local
  indices, and runs are split at tile column seams. Only non-empty tiles
  are returned, each with its offset and its own `GridSpec`. A property
  test checks every tile against `crop` record for record, and the
  reassembled tiles against the whole raster bit for bit
  (`tests/tiles.rs`). The README's Chunking recipe used to loop `crop`
  over the tiles, which scans every record for every tile; it now uses
  `slice_tiles`.
- **`crop` is public, and 0-based.** `src/crop.rs` (the Rust analogue of R's
  `crop_burn()` and Python's `BurnResult.crop()`) was never declared in
  `lib.rs`, so it was not compiled, and it still used the 1-based,
  inclusive indices of the tables before the switch to 0-based, half-open
  runs. It is now `controlledburn::crop`, re-bases to the window's top-left
  cell as row 0, col 0, and keeps runs half-open. A property test checks
  that cropping then materializing equals slicing the full raster, bit for
  bit, on random windows of a non-dyadic grid (`tests/crop.rs`). The book's
  crop chapter is now in its table of contents, so its examples are tested.
- **Canonical tables for the C++ core.** `tests/canon.rs` generates 240
  lattice-biased cases (rings with vertices on corners, walls and centre
  lines, holes, multipolygons, lines, points) on non-dyadic grids and keeps
  `fixtures/canon-tables.txt`, every record of both modes in parity mode,
  current. The C++ core checks itself against the same file
  (hypertidy/controlledburn `cpp/tests/test_canon_tables.cpp`): at this
  commit all 6737 records agree, fractions as bits. The file freezes the
  inputs as well as the outputs: the check re-burns the recorded grid and
  WKB bits rather than regenerating them, because the generator uses the
  platform's `sin` and `cos`, which differ in the last bit between Linux,
  macOS and Windows.

## 0.2.4

- **`Measure`: spherical area for grids whose rows are parallels.**
  `CoverageWindow::with_measure` / `set_measure` with
  `Measure::CosLatitude { phi0, dphi }` (latitude `phi0 + dphi * y`, radians)
  measures area with the element `cos(phi) dx dy`, as on the sphere for any
  regular longitude/latitude grid. Each edge piece books closed-form
  integrals of `cos(phi) dy` and `(x - k) cos(phi) dy` (x is linear in y
  along a piece), so the row accumulation is unchanged; coverage is the
  covered measure over the cell's (`Measure::cell_area`), and
  `Measure::ring_area` gives a ring's. `Measure::Planar` is the default and
  its arithmetic is unchanged. Tested against closed forms, quadrature and
  per-cell clipping (`tests/measure.rs`). For rwarp's conservative weights
  (rwarp `design/conservative.md`, stage 3). The C++ core does not have it
  yet.

## 0.2.3

- **Coverage: cells overcounted where an edge passes exactly through a
  cell corner, or where two edges nearly coincide** (#6). Found by the
  property suite; the C++ core had the same bugs (hypertidy/controlledburn#25).
  - An edge through a corner gave a crossing an ulp off the corner, and the
    next cell a sliver traversal along its wall where the ring only touches
    the corner. The walker now exits through the corner itself when the
    crossing is within the lattice tolerance of it (1e-9 of a cell).
  - A cell crossed by several traversals got its coverage from
    exactextract's chain chaser, which depends on the order of the
    traversal ends around the cell. Rounding can swap two nearly coincident
    ends (that sliver, or the two edges of a ring doubling back along
    itself), and the chaser then returns the complement: a whole cell too
    much. Such cells now use Green's theorem, the identity behind
    `CoverageWindow`: each traversal adds its own term, plus the cell
    width times the covered length of the left wall, read from the
    traversals to its left in the row. So does a single traversal lying
    along the cell walls, which does not split the cell.
  - A zero-width spike out to a corner vertex and back left a coverage of
    1e-16 instead of 0, whose record blanked a fully covered cell. A
    coverage below the sweep's 1e-6 tolerance now gets a record only when
    the ring passes through the cell.

  Before, three rings from #6 burned one cell each at 1.0 for 0.094, 0.739
  for 0.455 and 0.839 or 1.0 (depending on the start vertex) for 0.339; a
  sliver of area 4e-13 covered 7 to 11 whole cells at 2x resolution. On
  `main`, every 100,000-case property run failed (20 of 20). The golden
  table for `bowtie` gains six cells, along its closing edge on a column
  wall: it covered 26 cells, now the exact 32. Every other golden table
  is unchanged.

- **Approx: a cell-centre tie decided by rounding** (#6). Which side of a
  cell centre an edge passes is now an exact orientation test on the
  edge's endpoints, not a comparison of the rounded intercept with the
  centre. Two edges of a ring doubling back along a line both passed
  1e-13 to the right of a centre, one intercept rounded onto it, the
  left-inclusive rule counted it as left, and an empty cell burned.
  Intercepts are now ordered by column and side of the centre. The
  sliver above still burns two cells in approx mode: one centre is
  strictly inside it and the other is on its edge, both checked in
  rational arithmetic, so that is the centre rule's answer, not a tie.
  The property "a cell untouched in coverage is absent in approx" now
  exempts cells whose centre is within 1e-3 cells of an edge: coverage
  drops fractions below 1e-6, so there both can be right.

- **Coverage: a cell with edges only along its walls was filled.** With
  lattice-aligned input (polygons traced from a raster at half the grid's
  cell size, say), a cell can have polygon edges along both of its walls,
  one crossing the row centre up and the other down. It has zero coverage
  and zero net winding, so it got no boundary record, and the sweep filled
  it with the winding from its left: one whole cell too many. The cell now
  gets a boundary record whenever the ring passes through it, from one
  point on its boundary to another (not for corner touches or zero-width
  spikes, which enclose nothing). Found on DEA Waterbodies on a 50 m grid; the minimal case
  `POLYGON ((0 0, -2 5, 2 0, 1 0, 1 1, -1 3.5, 0 0.5, 0 0))` on unit cells
  covered 3.5 cells instead of 2.5. In parity mode, a zero-coverage cell
  inside the exterior on a hole's wall is now its own one-cell run, so
  the golden tables for `hole_ccw`, `hole_cw` and `nested_holes_multi`
  split one run each (cell values unchanged, coalesced runs unchanged).
  Same fix as the C++ core (hypertidy/controlledburn#29).

## 0.2.2

- **`CoverageWindow` and `RowPrefix`: exact coverage and exact integrals
  for small rings, without the traversal tables.** Both work in cell index
  space and split each ring edge at cell walls. `CoverageWindow` books
  each piece's area and height into a dense, reusable window and resolves
  coverage with one prefix pass per row (the cell accumulation of
  FreeType's "gray" rasteriser, kept exact). `RowPrefix` holds row prefix
  sums of one or more fields, so that a ring's weighted sum is a boundary
  integral: work proportional to the perimeter, and an edge shared by two
  rings integrated once. They agree with `burn` in coverage mode to its
  `f32` fractions (property-tested). Made for rwarp's exact footprint, one
  small ring per destination pixel, where `burn`'s per-call setup and
  per-cell records cost 1 to 2 microseconds a pixel; the window costs 0.2
  to 0.5. See `src/window.rs` for the identity and a C++ port note.

- **Coverage: holes removed too much area.** A cell crossed only by a
  hole's boundary, inside the exterior, got the hole's negative fraction
  alone and was dropped, instead of 1 minus the hole's fraction. A
  9 x 9 square with a 1.6 x 1.6 hole covered 77 cells instead of 78.44;
  on random valid polygons with holes, 98% had at least one wrong cell.
  The row sweep now tracks winding per ring and adds the cover of every
  ring that encloses a boundary cell without passing through it, kept
  incrementally so the cost doesn't grow with the number of holes.
  Polygons without holes are unchanged. The C++ core has the same bug
  (hypertidy/controlledburn).

- **Run coalescing dropped cells from multipolygons.** A multipolygon's
  runs were merged as if they were sorted by row and column, but they are
  only sorted within each part. When a later part had a run on a row to
  the left of an earlier part's last run on that row, the merge swallowed
  it, and those fully covered cells were lost (both modes). The runs are
  now sorted before merging, so the coalesced runs always equal the parity
  runs passed through `coalesce_runs`, and touching parts merge into one
  run. Parity mode (`coalesce_runs = false`) was not affected, nor is the
  C++ core, which does not coalesce. A multipolygon's runs now come out
  ordered by row and column instead of part by part. Found burning North
  Carolina counties onto a rotated grid (#7), where Craven county lost
  two cells.

  The book's one-degree world example was affected too: Canada lost five
  cells on row 16, so that burn now has 1,266 runs instead of 1,265, and
  `world-1deg.png` is regenerated.

## 0.2.1

Bug fixes only; no public API changes. Output changes only where the
previous output was wrong.

- **Coverage: vertices a few ulp off a cell edge** (#1). A vertex parsed
  from input (0.3) and a cell edge computed from the grid (3 * 0.1) can
  differ by a few ulp; the walk then followed the extrapolated segment off
  the grid and the polygon came back silently empty. Such vertices are now
  snapped onto the lattice before the walk, a crossing beyond a segment's
  end is clamped to the cell, and a ring walk that does not converge
  returns an error note instead of growing without bound. Same fix as
  hypertidy/controlledburn#16.
- **Coverage: a row-centre crossing on a column wall** (#1). Such a
  crossing registered no winding in either adjacent cell, and the row lost
  its interior runs from that column on. It now counts exactly once. Same
  fix as hypertidy/controlledburn#19.
- **Approx: a vertex on a row centre line** (#2). A vertex an ulp off the
  centre line as recomputed from the extent could flip the sign of its
  edge's intercept, and the row filled to the right edge of the grid. The
  crossing direction now comes from the edge's own endpoints.
- **Grid validation** (#3). A grid whose bounds are finite but whose span
  overflows `f64` (e.g. `-1e308` to `1e308`) gave an infinite cell size
  and an empty, `Ok` result. `GridSpec::validate()` now rejects it, as the
  C++ core does since hypertidy/controlledburn#23.
- **Tests.** A property-based suite (`tests/properties.rs`) over random
  rings with vertices steered onto cell corners, walls and centre lines;
  the shared fixtures and the golden C++ parity file gain the
  `wall_crossing` case (26 golden cases, 11 shared fixtures).
- **Docs.** `GridSpec` said row 1 is the top row; it is row 0. `NOTICE`
  now carries the MIT notice for the fasterize-derived approx scanline
  and pixel functions.

## 0.2.0

- **Breaking:** output tables are 0-based with half-open runs and `id = k`
  (0.1.0 was 1-based with inclusive `col_end` and `id = k + 1`).
- The book, the `geo-types` feature and the figure renderer.

## 0.1.0

- Initial pure-Rust port of the controlledburn C++ core.
