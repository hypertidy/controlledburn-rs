# DE-9IM per cell from the coverage traversal

Status: discussion. Scope: the proposal to decode DE-9IM(P, C) for a polygon P and each grid
cell box C from the coverage fraction f, a one-byte contact mask per edge
cell, and one bit per geometry; with a parallel `edge_contact: Vec<u8>`
column and a `touches` table (cells with f = 0 and wall or corner
contact). controlledburn-rs is the canonical implementation; C++ notes at
the end.

## Summary

- Related to controlledburn, not a separate topic: everything comes from
  the coverage walk, and nothing needs a second geometry.
- The plan is sound and the existing tables are untouched. Two
  corrections to the decoding table, one redundancy, and one real
  obstacle (the touches table cannot come from the current walk as is).
- Useful as: exact cell predicates (covers / intersects / touches /
  within per cell), an exact ALL_TOUCHED, and a sound prefilter for
  polygon-polygon predicates on a shared grid (measured below).
- Overbuild check: nobody else does cell-box DE-9IM from a traversal.
  geo (`Relate`), GEOS and DuckDB spatial do exact vector-vector DE-9IM;
  PostGIS raster has ST_Intersects(raster, geom) but via polygonising or
  pixel centres. The prefilter idea has prior art (raster-interval
  approximations for spatial joins, e.g. APRIL by Georgiadis, Mamoulis et
  al.), which used its own rasteriser; controlledburn's runs + edges are
  already its two lists (FULL, PARTIAL).

## The decoding table, checked

Rows are P (I, B, E), columns are C (I, B, E). P valid, C closed box.

| Entry | Proposal | Verdict |
|---|---|---|
| II | 2 if f > 0 else F | Correct (open sets meeting have positive area). |
| EI | 2 if f < 1 else F | Correct (P closed, full area implies P covers C). |
| EE | always 2 | Correct. |
| IE, BE | 2 / 1 unless P inside one cell | Correct; the bit must also say P is inside the grid, since a polygon with one in-grid cell may extend past the extent. |
| BI | mask bit | Redundant for valid P: BI = 1 iff 0 < f < 1. f = 1 means int C is in int P; f = 0 means a boundary through int C would leave area on one side. Keep the bit only if invalid input (zero-width spikes, slivers) should be reported, which is a validity flag, not topology. |
| BB | 1 / 0 / F from mask | Correct, and this is the one entry that genuinely needs new data. |
| IB, EB | mask bits | Mostly derivable: f = 0 gives IB = F; f = 1 gives EB = F. For 0 < f < 1 both are needed (P strictly inside C has IB = F). Values are only 1 or F. EB = F with f = 0 happens when a hole equals C, a nice edge case for tests. |

Run cells: the claim that all run cells share `212F01FF2` is wrong twice.
- A cell strictly inside P has BB = F, so the matrix is `212FF1FF2`.
- Run cells are not uniform. A boundary along a cell wall gives f = 1
  with BB = 1 (and IB = 1 on the other walls); P equal to C gives
  IB = F, BB = 1, IE = F. Checked: a 4x4 wall-aligned square burns to
  runs only, no edges.
- The predicate is the same (P covers C), so for "covers / contains"
  queries runs still need no storage. For the full matrix, wall-aligned
  run cells need a mask too. Cheapest fix: emit wall-aligned boundary
  cells as edge-table rows with f = 1 when contact output is on, or give
  the touches table a flag for "f = 1 side" so both sides of a wall
  appear.

## Where the data comes from

- The walker already knows entry and exit sides and the coordinates of
  each traversal, so IB / EB / BB for edge cells are a per-traversal
  bitmask computed in the same pass. Cost is one byte per edge record.
- The touches table is the obstacle. A segment lying on a shared wall
  is walked through one cell only; a vertex on a corner touches up to
  four cells, of which the walk visits one. The f = 0 neighbours are
  never visited. They have to be emitted explicitly when a traversal
  runs along a wall or passes through a corner, which is a small change
  in the walker but outside the exactextract-derived logic, so parity
  tests must check it does not alter `edges`.
- Snapping: the walker moves vertices within a relative tolerance of a
  wall onto the wall (checked: a square offset by 1e-14 burns
  identically to the aligned one). That is the right choice for
  coverage, but it means BB and the touches table are "contact within
  snapping tolerance". Document it; it is the same as GEOS with a
  precision model, not exact predicates.

## Non-affine grids

The walk is in index space, so the matrix is DE-9IM in index space.
- Affine and rotated: exact, the map is an affine bijection.
- Rectilinear: exact when edges are cut at the walls (the
  grid-frontends exact method). Vertex-only transforms bend edges in
  index space and can lose or invent corner contacts.
- Curvilinear: II / EI with f clearly away from 0 and 1 are fine; BB,
  IB, EB and touches are exactly the entries the 1-4 percent of a cell
  index-space error disturbs. Report them as index-space topology.
- The good case: polygons made from the grid itself (polygonised masks,
  cell outlines, coverage-aligned boundaries) lie on walls exactly in
  index space on any grid, and that is where BB and touches matter
  most.

## What it is for

1. Cell predicates as sparse masks: cells covered by P (runs), cells
   intersecting P (runs + edges + touches), cells touching only
   (touches), boundary cells (edges). Exact ALL_TOUCHED is "runs +
   edges + touches"; GDAL's all-touched misses wall and corner contacts
   and is a natural oracle comparison.
2. A prefilter for polygon-polygon predicates on a shared grid, from f
   alone (no mask needed). Experiment below.
3. Not exact polygon-polygon DE-9IM. Combining two cell matrices gives a
   three-way answer (true / false / refine); exact refinement belongs to
   geo or GEOS on the undecided pairs, ideally only on the shared edge
   cells.

## Experiment: prefilter soundness and yield

`examples/de9im_prefilter.rs` (`cargo run --release --example
de9im_prefilter`; geo is a dev-dependency). 300 random star polygons on a 100 x 100 domain, 6334 bbox-overlapping
ordered pairs, truth from `geo::Relate`. Rules use only runs (FULL) and
edges (PARTIAL):

- intersects true if a FULL cell of one meets any cell of the other;
  false if no cell of A is within one cell of a cell of B (the one-cell
  dilation is required because wall-aligned touching has no shared
  cell; with a touches table the dilation goes away).
- interiors meet true by the same rule; false if no shared cell.
- within (A in B) false if A has a FULL cell B does not cover fully, or
  a cell B does not have; true if every A cell is FULL in B.

| grid | intersects decided | interiors meet decided | within decided | wrong |
|---|---|---|---|---|
| 32 | 41.6% | 51.2% | 95.5% | 0 |
| 128 | 76.8% | 85.2% | 98.9% | 0 |
| 512 | 93.0% | 96.1% | 99.8% | 0 |
| 2048 | 98.1% | 99.0% | 99.9% | 0 |

The example runs 32 to 512; the 2048 row came from the same code with a
larger grid.

Burning all 300 polygons at 2048 took about 0.1 s. The filter in the
example uses hash sets and is slow; a real one merges sorted run
intervals per row, which is the APRIL design.

## Recommendation

Worth pursuing, in two small steps in controlledburn-rs:
1. `edge_contact: Vec<u8>` behind an option (IB, EB, BB bits; drop BI or
   keep it as a validity flag), plus wall-aligned f = 1 cells surfaced.
   Fixtures: wall-aligned square, corner-touching squares, hole equal to
   a cell, polygon inside one cell, polygon crossing the extent.
2. The `touches` table, emitted from wall-running and corner-passing
   traversals, with golden parity on all existing tables.
A predicates module (cell masks, ALL_TOUCHED, the pairwise prefilter) can
sit on top later, possibly in wxy rather than in the crate.

C++ port: the walker change is the same function in both (the C++
walk_polyline); keep the mask computation in its own helper so the port
is a copy, and extend the golden JSON with the two new columns so the
existing cross-language parity test covers them.
