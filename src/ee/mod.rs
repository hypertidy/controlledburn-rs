//! Geometry primitives derived from exactextract (Daniel Baston,
//! ISciences LLC, Apache-2.0): axis-aligned boxes with side and crossing
//! logic, perimeter-distance parameterisation, area measures and the
//! left-hand-area chain chaser. Ported with the arithmetic kept in the
//! same order so results are bit-identical to the C++.

pub mod bbox;
pub mod measures;
pub mod perimeter;
pub mod side;
pub mod traversal_areas;

