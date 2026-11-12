//! Minimal planar geometry model.
//!
//! The rasterizer needs only coordinates, ring structure and geometry kind:
//! no topology, no validity checking, no CRS. Rings are flat coordinate
//! vectors; orientation is derived from the signed area at burn time, so
//! input orientation does not matter.
//!
//! `GeometryCollection` is intentionally absent: mixed-dimension input
//! would produce a sparse table with inconsistent weight semantics (50 m^2
//! of polygon vs 50 m of line, indistinguishable in the output). Split
//! collections into homogeneous groups upstream.

/// A 2-D coordinate.
#[derive(Clone, Copy, Debug, PartialEq, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Coord {
    pub x: f64,
    pub y: f64,
}

impl Coord {
    #[inline]
    pub const fn new(x: f64, y: f64) -> Self {
        Coord { x, y }
    }

    #[inline]
    pub fn is_finite(&self) -> bool {
        self.x.is_finite() && self.y.is_finite()
    }
}

impl From<(f64, f64)> for Coord {
    #[inline]
    fn from((x, y): (f64, f64)) -> Self {
        Coord { x, y }
    }
}

impl From<[f64; 2]> for Coord {
    #[inline]
    fn from([x, y]: [f64; 2]) -> Self {
        Coord { x, y }
    }
}

/// An ordered coordinate sequence: a ring or a linestring.
pub type CoordSeq = Vec<Coord>;

/// One polygon: exterior ring first, then zero or more hole rings.
#[derive(Clone, Debug, PartialEq, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Polygon {
    pub rings: Vec<CoordSeq>,
}

impl Polygon {
    pub fn new(rings: Vec<CoordSeq>) -> Self {
        Polygon { rings }
    }

    pub fn exterior(&self) -> Option<&CoordSeq> {
        self.rings.first()
    }
}

/// A geometry of one of the six supported kinds.
///
/// Multi variants may hold zero or one part; the singular variants are
/// provided so the kind of the input is preserved (for example for
/// round-tripping), the burn treats both identically.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Geometry {
    Point(Coord),
    MultiPoint(Vec<Coord>),
    LineString(CoordSeq),
    MultiLineString(Vec<CoordSeq>),
    Polygon(Polygon),
    MultiPolygon(Vec<Polygon>),
}

impl Geometry {
    /// True when the geometry holds no coordinates at all. A point with
    /// both coordinates NaN is WKB's `POINT EMPTY`, and counts as none.
    pub fn is_empty(&self) -> bool {
        match self {
            Geometry::Point(c) => is_empty_point(c),
            Geometry::MultiPoint(p) => p.iter().all(is_empty_point),
            Geometry::LineString(l) => l.is_empty(),
            Geometry::MultiLineString(l) => l.is_empty(),
            Geometry::Polygon(p) => p.rings.is_empty(),
            Geometry::MultiPolygon(p) => p.is_empty(),
        }
    }

    /// The point parts (one for `Point`, all for `MultiPoint`, none otherwise).
    pub fn points(&self) -> &[Coord] {
        match self {
            Geometry::Point(c) => std::slice::from_ref(c),
            Geometry::MultiPoint(p) => p,
            _ => &[],
        }
    }

    /// The line parts.
    pub fn lines(&self) -> &[CoordSeq] {
        match self {
            Geometry::LineString(l) => std::slice::from_ref(l),
            Geometry::MultiLineString(l) => l,
            _ => &[],
        }
    }

    /// The polygon parts.
    pub fn polygons(&self) -> &[Polygon] {
        match self {
            Geometry::Polygon(p) => std::slice::from_ref(p),
            Geometry::MultiPolygon(p) => p,
            _ => &[],
        }
    }

    /// Iterate every coordinate of the geometry.
    pub fn coords(&self) -> impl Iterator<Item = &Coord> {
        let pts = self.points().iter();
        let lines = self.lines().iter().flatten();
        let polys = self.polygons().iter().flat_map(|p| p.rings.iter().flatten());
        pts.chain(lines).chain(polys)
    }

    /// True when every coordinate is finite, empty points (both
    /// coordinates NaN) aside.
    pub fn is_finite(&self) -> bool {
        let pts = self.points().iter().all(|c| c.is_finite() || is_empty_point(c));
        let lines = self.lines().iter().flatten().all(Coord::is_finite);
        pts && lines && self.polygons().iter().flat_map(|p| p.rings.iter().flatten()).all(Coord::is_finite)
    }
}

/// WKB's `POINT EMPTY`: both coordinates NaN.
fn is_empty_point(c: &Coord) -> bool {
    c.x.is_nan() && c.y.is_nan()
}

/// Signed area of a closed ring (shoelace). Positive = counter-clockwise
/// in a conventional (x east, y north) frame.
///
/// The summation is anchored at the first vertex, in the same order as
/// the C++ core, so results are bit-identical.
pub fn signed_area(ring: &[Coord]) -> f64 {
    if ring.len() < 3 {
        return 0.0;
    }
    let mut sum = 0.0;
    let x0 = ring[0].x;
    let mut i = 1;
    while i + 1 < ring.len() {
        let x = ring[i].x - x0;
        let y_next = ring[i + 1].y;
        let y_prev = ring[i - 1].y;
        sum += x * (y_next - y_prev);
        i += 1;
    }
    sum / 2.0
}

/// True when the ring is counter-clockwise (positive signed area).
#[inline]
pub fn is_ccw(ring: &[Coord]) -> bool {
    signed_area(ring) > 0.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shoelace_sign() {
        let ccw =
            vec![Coord::new(0., 0.), Coord::new(1., 0.), Coord::new(1., 1.), Coord::new(0., 1.), Coord::new(0., 0.)];
        assert_eq!(signed_area(&ccw), 1.0);
        let mut cw = ccw.clone();
        cw.reverse();
        assert_eq!(signed_area(&cw), -1.0);
        assert!(is_ccw(&ccw) && !is_ccw(&cw));
        assert_eq!(signed_area(&ccw[..2]), 0.0);
    }
}

/// A geometry that cannot be represented: a collection (from the
/// `geo-types` or `geo-traits` conversions).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeometryCollectionError;

impl std::fmt::Display for GeometryCollectionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GeometryCollection is not supported; split into homogeneous groups")
    }
}

impl std::error::Error for GeometryCollectionError {}
