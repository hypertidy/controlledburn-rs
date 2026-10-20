# Topology and meshes: what the burn already knows

Status: discussion. Follows `de9im-cells.md`. Relates the coverage
traversal to mesh representations (shared vertices and edges, as in
arc-node topology, triangulations and silicate-style vertex and edge
tables), and sets out what controlledburn could expose or accept if
topology is made explicit.

## Summary

- A regular grid is a quad mesh. Burning a polygon onto it is the
  overlay of two meshes: the polygon's boundary and the lattice of cell
  walls. The walker already computes that overlay, piece by piece: it
  cuts each ring at cell walls and records, for each piece, the cell it
  lies in and the sides it enters and leaves by.
- The per-cell DE-9IM tables (`edge_contact`, `touches`) are incidence
  tables of that overlay. Stored as an overlay (pieces with a cell id
  and a left and right polygon id), they become views rather than
  separate products.
- For features that already share a mesh, spatial predicates between
  them need no geometry: shared edge means touch along a line, shared
  vertex only means touch at a point, and a partition has disjoint
  interiors by construction.
- Explicit topology does not remove the hard part, which is building it
  from independent simple-features data (noding, snapping). That stays
  a conversion step.
- Prospects for controlledburn-rs, in order of value: (1) a shared-arc
  input mode for coverages and meshes, (2) an overlay output (boundary
  pieces per cell with left and right ids), (3) the DE-9IM cell tables
  as views of (2). Each is opt-in and leaves the existing tables
  unchanged.

## The burn as a mesh overlay

Three meshes are in play.

| mesh | vertices | edges | faces |
|---|---|---|---|
| grid | cell corners | cell walls | cells |
| features | ring vertices | ring segments (shared, if a coverage) | polygons |
| overlay | both, plus wall crossings | segments cut at walls, wall pieces | polygon-in-cell pieces |

The walker produces the overlay's boundary edges (traversals: a piece of
ring inside one cell, with entry and exit sides). The left-hand-area
chain chaser closes each cell's pieces against the cell perimeter, which
is building the overlay faces and measuring them. The winding sweep
labels the cells the boundary never enters (runs). Nothing in the
overlay is new; only the incidence is discarded after the areas are
taken.

That incidence is what the DE-9IM cell matrix needs:

- BB (boundary meets cell boundary): a piece starts, ends or runs on a
  wall or corner.
- IB, EB: which arcs of the cell perimeter lie between pieces on the
  inside or outside, which the chain chaser walks anyway.
- `touches`: faces of the overlay with zero area on one side of a wall
  piece, that is the cell across the wall from a piece lying on it.

## Features that already share a mesh

When features come as one mesh (a polygon coverage with shared arcs, a
triangulation, a silicate-style model), relations between features are
table lookups:

| relation | from the mesh |
|---|---|
| touches along a line | the two faces share an edge |
| touches at a point | they share a vertex but no edge |
| disjoint | neither |
| interiors overlap | impossible in a partition; for overlapping layers, only via an overlay |
| within, contains | face nesting, if the mesh records it (holes, islands) |

This is exact and needs no snapping tolerance, because sharing is by
identity rather than by coordinate equality. It is what a topology
model (arc-node, PostGIS topology, TopoJSON) gives, and why the
simple-features predicates in GEOS, geo or DuckDB spatial have to
rediscover it with noding and robust arithmetic.

Feature-to-cell relations are still geometry, because the grid's walls
are not edges of the feature mesh until the overlay is computed. That
is the part controlledburn does.

## Prospect 1: shared-arc input

Burning a coverage today means burning each polygon's rings, so every
interior arc is walked twice, once in each direction. A shared-arc mode
takes arcs with a left and a right polygon id:

- Walk each arc once. Each traversal is credited to the left polygon as
  is, and to the right polygon reversed.
- Per cell, per polygon, the chain chaser runs on that polygon's pieces,
  exactly as now. Winding deltas go to both sides with opposite signs.
- Holes, islands and the outside are ids like any other (the outside is
  a reserved id that is not emitted).

What it gives:

- Up to half the walking for a coverage (every interior arc once instead
  of twice). The chain chasing per polygon per cell is unchanged.
- One set of crossing coordinates per arc for both neighbours, so
  adjacency and the `touches` table are consistent by construction
  rather than by two walks agreeing.
- A natural input for meshes: a triangulation's edge table with left and
  right triangle ids is exactly this.

What it does not give: a fix for gaps or overlaps between neighbours.
`examples/partition_sum.rs` burns a jittered 24 x 24 quad partition
polygon by polygon and sums the fractions in each cell:

| grid | edge records | max abs(sum - 1) |
|---|---|---|
| 50 x 50 | 5477 | 9.3e-7 |
| 173 x 173 | 19775 | 4.1e-8 |
| 512 x 512 | 58995 | 1.0e-6 |

Every cell is reached, and the sums are 1 to the single precision the
fractions are stored in. Burning shared edges twice is already
consistent in area; the gain of the arc mode is topology and work, not
area.

Cost: a new input type (arcs plus a face table, or a mesh with left and
right ids), a validity check that the arcs close into rings per id, and
a parity test that the arc mode reproduces the per-polygon tables on the
existing fixtures.

## Prospect 2: overlay output

Optionally return the overlay pieces themselves: for each boundary cell,
the ring pieces inside it with their cell id and left and right ids
(coordinates as an arena plus offsets, as the walker holds them now).

What it gives:

- Exact polygon-in-cell geometry, not just its area. That is the input
  for area-weighted interpolation and conservative regridding between a
  grid and polygons or a mesh, the same weights ESMF computes with its
  mesh-to-mesh overlay.
- Polygonising back: cell runs plus pieces reconstruct the polygon on the
  grid, which is a round-trip test.
- The DE-9IM cell tables computed from the pieces, rather than as a
  second set of flags in the walker.

Cost: memory proportional to the boundary, which is the same order as
the existing `edges` table, and an output type to version.

## Prospect 3: DE-9IM cell tables as views

With prospect 2, `edge_contact` and `touches` are derived, and their
rules can live in one place:

- `edge_contact` bits from where pieces meet the cell perimeter.
- `touches` from pieces lying on a wall or corner, emitted for the cell
  on the other side.
- Run cells with a wall-aligned boundary get their BB = 1 from the same
  pieces, which closes the gap noted in `de9im-cells.md`.

Without prospect 2 the walker computes the same bits directly, which is
cheaper if nothing else needs the pieces. Either is fine; the point is
that both read the same incidence.

## Non-regular grids

The overlay is computed in index space, so on non-affine grids it is the
overlay with the cell walls straight in index space. For a curvilinear
grid in world space the cells are bilinear quads with bent sides, and
the world-space clip in `grid-frontends.md` is the same overlay with the
other walls. Seen as meshes, the two methods differ only in which mesh
supplies the cell walls; both produce pieces with cell and polygon ids.
For meshes built from the grid itself (polygonised masks, cell
outlines), the shared edges are cell walls and the overlay is exact in
index space on any grid.

## What stays outside controlledburn

- Building topology from independent features: noding, snapping,
  repairing. That belongs to conversion tools (GEOS, geo, or a dedicated
  topology builder), with controlledburn taking the result.
- Predicates between features of one mesh: lookups on the mesh, not a
  burn.
- Feature-feature DE-9IM for features from different meshes: the cell
  prefilter in `de9im-cells.md` decides most pairs; the rest go to an
  exact engine.

## C++ port

The walker, the chain chaser and the winding sweep are the same
functions in both implementations. The arc mode changes the driver
(which rings feed the walk and how traversals are credited), not those
functions, so the port is a new driver plus the input type. The overlay
output is a copy of what `walk_polyline` already holds. Extend the golden
JSON with an arc-mode case built from an existing fixture, so the
cross-language parity test covers both.
