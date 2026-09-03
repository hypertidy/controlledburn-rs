// Derived from exactextract measures.cpp
// Copyright (c) 2018 ISciences, LLC. Apache License 2.0.

use crate::geometry::Coord;

/// Signed area with exactextract's sign convention (clockwise positive),
/// the OPPOSITE of `crate::signed_area`. Only its absolute value is used
/// by the coverage math; it is kept separate so the summation order
/// matches the C++ exactly.
pub fn area_signed(ring: &[Coord]) -> f64 {
    if ring.len() < 3 {
        return 0.0;
    }
    let mut sum = 0.0;
    let x0 = ring[0].x;
    let mut i = 1;
    while i < ring.len() - 1 {
        let x = ring[i].x - x0;
        let y1 = ring[i + 1].y;
        let y2 = ring[i - 1].y;
        sum += x * (y2 - y1);
        i += 1;
    }
    sum / 2.0
}

#[inline]
pub fn area(ring: &[Coord]) -> f64 {
    area_signed(ring).abs()
}

/// Sum of segment lengths.
#[allow(dead_code)]
pub fn length(coords: &[Coord]) -> f64 {
    let mut sum = 0.0;
    for w in coords.windows(2) {
        let dx = w[1].x - w[0].x;
        let dy = w[1].y - w[0].y;
        sum += (dx * dx + dy * dy).sqrt();
    }
    sum
}
