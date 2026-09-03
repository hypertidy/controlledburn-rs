# controlledburn

Sparse scanline rasterization of polygons, lines and points onto a
regular grid, with exact coverage fractions for polygon boundary cells.
O(perimeter) in time and memory: no pixel buffer is ever materialized by
the burn itself.

Pure Rust port of the C++ core of
[controlledburn](https://github.com/hypertidy/controlledburn). No
dependencies, no unsafe code, output bit-identical to the C++ (see
Parity below).

```rust
use controlledburn::{burn, BurnOptions, Coord, Geometry, GridSpec, Polygon};

let square = Geometry::Polygon(Polygon::new(vec![vec![
    Coord::new(2.5, 4.5), Coord::new(6.5, 4.5), Coord::new(6.5, 8.5),
    Coord::new(2.5, 8.5), Coord::new(2.5, 4.5),
]]));
let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10); // extent, ncol, nrow
let r = burn(&[square], &grid, BurnOptions::default()).unwrap();

// r.runs   interior cells:      (row, col_start, col_end, id)
// r.edges  boundary cells:      (row, col, fraction, id), fraction in (0, 1)
// r.lines  line cells:          (row, col, length, id), CRS units
// r.points point cells:         (row, col, id)
// r.notes  non-fatal problems per input geometry
assert!((r.covered_cells() * grid.dx() * grid.dy() - 16.0).abs() < 1e-5);
```

Indices are 1-based, row 1 is the top row, geometry `k` (0-based input
position) gets id `k + 1`. WKB input (ISO and EWKB, both byte orders,
Z/M skipped) goes through `burn_wkb`.

## Modes

- `BurnMode::Coverage` (default): exact coverage fractions via analytical
  cell traversal (derived from exactextract). Total burned coverage equals
  the exact polygon area inside the grid.
- `BurnMode::Approx`: the cell-centre rule (fasterize semantics), runs
  only, left-inclusive on ties. Much faster.

Lines and points are unaffected by mode. `materialize` reduces a result
into a dense `f64` buffer with fasterize-style pixel functions.

## Differences from the C++ core

- Unclosed polygon rings are closed by appending the first coordinate.
- Geometries with non-finite coordinates are skipped with a note (the
  C++ walker does not terminate on a NaN vertex).
- Adjacent runs on a row are coalesced by default;
  `BurnOptions::parity(mode)` turns that off for record-for-record
  equality with the C++ tables.
- Sorts are stable, so ties among equal-x intercepts in Approx mode are
  deterministic.

## Parity

`tests/golden.rs` compares every output table, in order, against
`fixtures/controlledburn-golden.json`, which was dumped from the C++
build for 25 cases (the 10 shared cross-language fixtures plus 15
edge-case probes) in both modes. Edge fractions and line lengths are
compared as f32 bit patterns. `tests/fixtures.rs` runs the shared
aggregate contract from `fixtures/{geometries,expected}.csv`.

## Performance

200,000-vertex jittered ring, single core (`cargo bench`):

| grid | mode | C++ core | this crate |
|---|---|---|---|
| 4096 x 4096 | Coverage | 3.23 s | 1.12 s |
| 4096 x 4096 | Approx | 69 ms | 77 ms |
| 1024 x 1024 | Coverage | 0.65 s | 0.23 s |

The Coverage-mode gain comes from data structures only (a hash-indexed
cell list instead of `std::map`, push-and-merge instead of a per-row
linear scan, a coordinate arena instead of one vector per traversal, and
an O(k log k) chain search in the left-hand-area routine with the same
tie-breaking); the arithmetic is unchanged.

## License

Apache-2.0. See NOTICE for the exactextract attribution.
