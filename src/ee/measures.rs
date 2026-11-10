// Derived from exactextract measures.cpp
// Copyright (c) 2018 ISciences, LLC. Apache License 2.0.

use crate::geometry::Coord;

/// Signed area with exactextract's sign convention (clockwise positive),
/// the OPPOSITE of `crate::signed_area`. Only its absolute value is used
/// by the coverage math; it is kept separate so the summation order
/// matches the C++ exactly.
pub fn area_signed(ring: &[Coord]) -> f64 {
    area_signed_by(ring.len(), |i| ring[i])
}

/// `area_signed` over a ring given by its length and an accessor, for a
/// ring assembled from pieces without collecting it (the same operations
/// in the same order, so the same bits).
#[inline]
pub fn area_signed_by(n: usize, at: impl Fn(usize) -> Coord) -> f64 {
    if n < 3 {
        return 0.0;
    }
    let mut sum = 0.0;
    let x0 = at(0).x;
    let mut i = 1;
    while i < n - 1 {
        let x = at(i).x - x0;
        let y1 = at(i + 1).y;
        let y2 = at(i - 1).y;
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
