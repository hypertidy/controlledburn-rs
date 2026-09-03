# Wrapping the crate

*Coming in a later revision.* What a binding needs: the four tables are
`#[repr(C)]` structs of `i32` and `f32`, so they cross an FFI boundary
as flat arrays; `burn_wkb` accepts the bytes a database or geometry
library already has; the `serde` feature serialises everything; the
`geo-types` feature converts from the Rust geospatial ecosystem's
common types. Sketches for pyo3 and extendr will live here.
