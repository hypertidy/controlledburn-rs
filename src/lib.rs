//! # controlledburn
//!
//! Sparse scanline rasterization of polygons, lines and points onto a
//! regular grid, with exact coverage fractions for polygon boundary
//! cells. O(perimeter) in time and memory: no pixel buffer is ever
//! materialized by the burn itself.
//!
//! This crate is a pure Rust port of the C++ core of
//! [controlledburn](https://github.com/hypertidy/controlledburn). It has
//! no dependencies and no unsafe code.
//!
//! ## The contract
//!
//! Input: geometries ([`Geometry`]), a grid ([`GridSpec`]: extent plus
//! column/row counts, row 1 at the top) and [`BurnOptions`].
//!
//! Output ([`BurnResult`]), four typed tables with 0-based indices (row 0
//! at the top) and half-open run ranges. Geometry `k` in the input slice
//! has `id = k`:
//!
//! | table    | columns                          | meaning                                  |
//! |----------|----------------------------------|------------------------------------------|
//! | `runs`   | row, col_start, col_end, id      | fully covered cells `col_start..col_end` |
//! | `edges`  | row, col, fraction (f32), id     | polygon boundary cells, fraction in (0,1)|
//! | `lines`  | row, col, length (f32), id       | line length inside the cell, CRS units   |
//! | `points` | row, col, id                     | one record per point inside the grid     |
//! | `notes`  | geom_index, message              | non-fatal problems, per input geometry   |
//!
//! Invariant (Coverage mode): `cell_area * (run cells + sum of edge
//! fractions)` equals the exact polygon area inside the grid.
//!
//! ## Modes
//!
//! * [`BurnMode::Coverage`] (default): exact fractions via analytical
//!   cell traversal, derived from exactextract.
//! * [`BurnMode::Approx`]: cell-centre rule (fasterize semantics), runs
//!   only, much faster.
//!
//! Lines and points are unaffected by mode.
//!
//! ## Example
//!
//! ```
//! use controlledburn::{burn, BurnOptions, Coord, Geometry, GridSpec, Polygon};
//!
//! let square = Geometry::Polygon(Polygon::new(vec![vec![
//!     Coord::new(2.5, 4.5), Coord::new(6.5, 4.5), Coord::new(6.5, 8.5),
//!     Coord::new(2.5, 8.5), Coord::new(2.5, 4.5),
//! ]]));
//! let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
//! let r = burn(&[square], &grid, BurnOptions::default()).unwrap();
//! let area = r.covered_cells() * grid.dx() * grid.dy();
//! assert!((area - 16.0).abs() < 1e-5);
//! ```
//!
//! ## Differences from the C++ core
//!
//! * Unclosed polygon rings are closed by appending the first coordinate.
//! * Geometries with non-finite coordinates are skipped with a note.
//! * Adjacent runs are coalesced by default; `BurnOptions::parity` turns
//!   that off for record-for-record equality with the C++ tables.
//! * Sorts are stable, so ties among equal-x intercepts in Approx mode
//!   are deterministic.
//! * Indices are 0-based with half-open run ranges and `id = k`; the C++
//!   core (and its R and Python bindings) emit 1-based indices, inclusive
//!   `col_end` and `id = k + 1`. `BurnOptions::parity` does not change
//!   this; the golden tests apply the offset when comparing.

#![forbid(unsafe_code)]
#![warn(missing_debug_implementations)]

mod approx;
mod burn;
mod coalesce;
mod coverage;
mod ee;
mod error;
#[cfg(feature = "geo-types")]
pub mod geo;
mod geometry;
mod grid;
mod line;
mod materialize;
mod output;
mod point;
mod polygon;
mod walker;
mod wkb;

pub use burn::{burn, burn_wkb, BurnMode, BurnOptions};
pub use coalesce::coalesce_runs;
pub use error::{BurnError, WkbError};
pub use geometry::{is_ccw, signed_area, Coord, CoordSeq, Geometry, Polygon};
pub use grid::GridSpec;
pub use materialize::{materialize, EdgePolicy, MaterializeOptions, PixelFn};
pub use output::{BurnResult, GridEdge, GridLine, GridPoint, GridRun, Note};
pub use wkb::parse_wkb;
