# Changelog

## Unreleased

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

## Unreleased

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
