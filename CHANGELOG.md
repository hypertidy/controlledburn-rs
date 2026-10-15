# Changelog

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
