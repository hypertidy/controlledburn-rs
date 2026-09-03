// Derived from exactextract perimeter_distance.cpp
// Copyright (c) 2018 ISciences, LLC. Apache License 2.0.

use super::bbox::BBox;
use crate::geometry::Coord;

/// Distance along the box perimeter, counter-clockwise from the
/// bottom-left corner: up the left side, along the top, down the right,
/// back along the bottom. `c` must lie on the boundary (the C++ throws
/// otherwise; here it is a debug assertion and the bottom-side formula is
/// used as a fallback).
pub fn perimeter_distance(b: &BBox, c: &Coord) -> f64 {
    let (x, y) = (c.x, c.y);
    if x == b.xmin {
        // Left
        return y - b.ymin;
    }
    if y == b.ymax {
        // Top
        return (b.ymax - b.ymin) + x - b.xmin;
    }
    if x == b.xmax {
        // Right
        return (b.xmax - b.xmin) + (b.ymax - b.ymin) + b.ymax - y;
    }
    debug_assert!(y == b.ymin, "perimeter_distance: coordinate not on boundary");
    // Bottom
    (b.xmax - b.xmin) + 2.0 * (b.ymax - b.ymin) + (b.xmax - x)
}

/// Counter-clockwise perimeter distance from `measure1` back to `measure2`.
/// Retained as the reference definition for the sorted chain search.
#[allow(dead_code)]
#[inline]
pub fn perimeter_distance_ccw(measure1: f64, measure2: f64, perimeter: f64) -> f64 {
    if measure2 <= measure1 {
        measure1 - measure2
    } else {
        perimeter + measure1 - measure2
    }
}
