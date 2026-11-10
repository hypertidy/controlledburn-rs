# controlledburn

Sparse scanline rasterization of polygons, lines and points onto a
regular grid, with exact coverage fractions for polygon boundary cells.
O(perimeter) in time and memory: no pixel buffer is ever materialized by
the burn itself.

Pure Rust port of the C++ core of
[controlledburn](https://github.com/hypertidy/controlledburn). No
dependencies, no unsafe code, output bit-identical to the C++ (see
Parity below).

- API reference: [docs.rs/controlledburn](https://docs.rs/controlledburn)
- The book, with real data: [hypertidy.github.io/controlledburn-rs](https://hypertidy.github.io/controlledburn-rs/)

```rust
use controlledburn::{burn, BurnOptions, Coord, Geometry, GridSpec, Polygon};

let square = Geometry::Polygon(Polygon::new(vec![vec![
    Coord::new(2.5, 4.5), Coord::new(6.5, 4.5), Coord::new(6.5, 8.5),
    Coord::new(2.5, 8.5), Coord::new(2.5, 4.5),
]]));
let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10); // extent, ncol, nrow
let r = burn(&[square], &grid, BurnOptions::default()).unwrap();

// r.runs   interior cells:      (row, col_start, col_end, id), col_start..col_end
// r.edges  boundary cells:      (row, col, fraction, id), fraction in (0, 1)
// r.lines  line cells:          (row, col, length, id), CRS units
// r.points point cells:         (row, col, id)
// r.notes  non-fatal problems per input geometry
assert!((r.covered_cells() * grid.dx() * grid.dy() - 16.0).abs() < 1e-5);
```

Indices are 0-based with row 0 at the top, run ranges are half-open,
and geometry `k` in the input has `id = k` (`id_offset + k` with
`BurnOptions::with_id_offset`). A multi-part geometry has one record per
cell: parts crossing the same cell are summed. WKB input (ISO and EWKB,
both byte orders, Z/M skipped) goes through `burn_wkb`.

## Modes

- `BurnMode::Coverage` (default): exact coverage fractions via analytical
  cell traversal (derived from exactextract). Total burned coverage equals
  the exact polygon area inside the grid.
- `BurnMode::Approx`: the cell-centre rule (fasterize semantics), runs
  only, left-inclusive on ties. Much faster.

Lines and points are unaffected by mode. `materialize` reduces a result
into a dense `f64` buffer with fasterize-style pixel functions.

[playproj: rasterize](https://hypertidy.github.io/playproj/rasterize.html) compares the
two rules interactively in the browser, on world countries and
coastlines in several projections, along with what each costs to store
in tiles; its "runs of whole cells plus a list of edge fractions" is
this crate's runs and edges tables.

## Chunking

A grid is four numbers and two counts: the extent `xmin, ymin, xmax,
ymax` and `ncol, nrow`. That fixes every cell edge, so a tile, chunk or
window is a range of rows and columns, and its extent is integer
arithmetic on the grid. The burn output is sparse and indexed on the
whole grid, so it is already chunkable: burn once, then `slice_tiles`
routes every record to its tile in one pass and you `materialize` each
tile. Tiles never need a second burn, and each geometry is walked once
however many tiles it spans.

```rust
let r = burn(&geoms, &grid, BurnOptions::default())?;
for t in slice_tiles(&r, &grid, 256, 256) { // tile shape in cells; empty tiles are skipped
    let mut buf = vec![f64::NAN; (t.grid.ncol * t.grid.nrow) as usize];
    materialize(&t.result, &mut buf, t.grid.ncol, t.grid.nrow, None, &opts)?;
    // write buf at cell offset (t.row_off, t.col_off), tile (t.tile_row, t.tile_col)
}
```

The tiles put back together are bit-identical to materializing the
whole grid (`tests/tiles.rs`). For a single window, or tile extents that
are not a regular grid, use `crop`. It returns the same tables but scans
the whole result on every call, so don't loop it over a tile grid.

- *Parallel*: split the geometries across workers, burn each subset
  onto the same grid with `BurnOptions::default().with_id_offset(start)`,
  and concatenate the tables. Records are per geometry, so nothing needs
  merging.
- *Output too large to hold*: burn a batch of geometries, slice it into
  tiles, merge each tile into the output, drop the batch. A tile can
  hold records from more than one batch, so merge rather than overwrite.
  Memory is bounded by the batch's perimeter, not by the grid.
- *Dense windows*: `CoverageWindow` resolves coverage straight into a
  window given by its integer row and column offsets, and a sub-window
  equals those cells of the whole (`tests/window.rs`).

## Features

- `serde`: `Serialize`/`Deserialize` on the public types.
- `geo-types`: `From`/`TryFrom` conversions from `geo_types` geometries.

## Differences from the C++ core

- Unclosed polygon rings are closed by appending the first coordinate.
- Geometries with non-finite coordinates are skipped with a note (the
  C++ walker does not terminate on a NaN vertex).
- Adjacent runs on a row are coalesced by default, and a multi-part
  geometry has one record per cell (the parts' fractions or lengths
  summed, a sum reaching a full cell becoming a run; overlapping parts
  are capped at a full cell, not unioned). `BurnOptions::parity(mode)`
  turns both off for record-for-record equality with the C++ tables.
- Rings of fewer than 4 coordinates once closed, and lines of one
  coordinate, are skipped with a note; the C++ core skips them silently.
  An all-NaN point (WKB `POINT EMPTY`) is empty, not non-finite.
- Indices are 0-based, run ranges half-open and `id = k`, the same as the
  C++ core and the Python package; the R package adds 1 (1-based indices,
  inclusive `col_end`, `id = k + 1`). The golden file predates the C++
  core's switch to 0-based and is in the 1-based form, so the golden
  tests apply that offset when comparing.
- Sorts are stable, so ties among equal-x intercepts in Approx mode are
  deterministic.

## Parity

`tests/golden.rs` compares every output table, in order, against
`fixtures/controlledburn-golden.json`, which was dumped from the C++
build for 31 cases (the 14 shared cross-language fixtures plus 17
edge-case probes) in both modes. Edge fractions and line lengths are
compared as f32 bit patterns. `tests/fixtures.rs` runs the shared
aggregate contract from `fixtures/{geometries,expected}.csv`.

## Properties

`tests/properties.rs` generates random simple rings on non-dyadic grids,
with vertices deliberately steered onto cell corners, cell walls, row and
column centre lines and the grid boundary (the places where parsed input
and computed cell edges disagree by an ulp), and checks invariants that
must hold for every input: coverage equals the ring's planar area clipped
to the grid; the tables are well formed; the cover matrix is invariant
under ring orientation and start vertex; it matches an exact per-cell
clipping oracle on small grids; approx and coverage agree on full and
empty cells; burning at 2x resolution and averaging 2x2 blocks equals
burning at 1x. CI runs 10000 cases per property. For a longer campaign:

```
PROPTEST_CASES=100000 cargo test --release --test properties
```

Shrunken failures are persisted in `tests/properties.proptest-regressions`
and replayed first on every run; commit that file when it changes. Rings
that touch themselves (a repeated vertex, a vertex on another edge) are
invalid polygons; the generator filters them and the crate gives no
defined result for them.

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
