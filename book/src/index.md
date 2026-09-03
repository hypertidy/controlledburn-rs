# controlledburn

`controlledburn` rasterizes polygons, lines and points onto a regular grid
and returns the result as four small tables instead of a pixel buffer.
For polygons it records which cells are fully inside (as runs) and, for
each boundary cell, the exact fraction of the cell the polygon covers.
Time and memory are proportional to the polygon perimeter, not to the
number of grid cells, so a grid with billions of cells is not a problem.

The crate is pure Rust with no dependencies and no unsafe code, and its
output is bit-identical to the C++ core it was ported from. It is meant
to be wrapped: the tables are plain `repr(C)` structs, WKB goes in
directly, and nothing in the API assumes a particular geometry library,
raster format or language on the other side.

This book is about the representation. Each chapter takes real data
(Natural Earth 110m countries, 177 features, public domain), burns it,
and does something with the tables that a dense raster would make slow,
large or impossible. Every code block on these pages is compiled and run
against the crate in CI; the numbers in the prose are asserted in the
code.

| chapter | what it shows |
|---|---|
| [The contract](contract.md) | one polygon, one grid, the four tables, the invariant |
| [Sparse versus dense](sparse-vs-dense.md) | the world at 1, 0.1 and 0.01 degrees: table sizes against dense array sizes |
| [Coverage versus approx](coverage-vs-approx.md) | exact fractions against the cell-centre rule |
| [Lines and points](lines-and-points.md) | the same tables for lower-dimensional input |
| [Materialize, chunk by chunk](materialize.md) | turning the tables into pixels one window at a time |
| [Zonal statistics from the tables](zonal.md) | area-weighted means straight from runs and edges |
| [Wrapping the crate](wrapping.md) | what a Python, R or database binding needs to know |

API reference: [docs.rs/controlledburn](https://docs.rs/controlledburn).
Source: [github.com/hypertidy/controlledburn-rs](https://github.com/hypertidy/controlledburn-rs).
