//! Conversions from `geo-types` (feature `geo-types`).
//!
//! `Point`, `MultiPoint`, `Line`, `LineString`, `MultiLineString`,
//! `Polygon`, `MultiPolygon`, `Rect` and `Triangle` convert with `From`.
//! `Geometry` converts with `TryFrom`, failing on `GeometryCollection`
//! (mixed dimensions have no place in the four-table output; split it
//! upstream).

use crate::geometry::{Coord, CoordSeq, Geometry, Polygon};

/// A `geo_types::Geometry` that cannot be represented: a collection.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeometryCollectionError;

impl std::fmt::Display for GeometryCollectionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GeometryCollection is not supported; split into homogeneous groups")
    }
}

impl std::error::Error for GeometryCollectionError {}

#[inline]
fn coord(c: geo_types::Coord<f64>) -> Coord {
    Coord::new(c.x, c.y)
}

fn line_string(l: &geo_types::LineString<f64>) -> CoordSeq {
    l.0.iter().map(|c| coord(*c)).collect()
}

fn polygon(p: &geo_types::Polygon<f64>) -> Polygon {
    let mut rings = Vec::with_capacity(1 + p.interiors().len());
    rings.push(line_string(p.exterior()));
    rings.extend(p.interiors().iter().map(line_string));
    Polygon { rings }
}

impl From<&geo_types::Point<f64>> for Geometry {
    fn from(p: &geo_types::Point<f64>) -> Self {
        Geometry::Point(coord(p.0))
    }
}

impl From<&geo_types::MultiPoint<f64>> for Geometry {
    fn from(p: &geo_types::MultiPoint<f64>) -> Self {
        Geometry::MultiPoint(p.0.iter().map(|q| coord(q.0)).collect())
    }
}

impl From<&geo_types::Line<f64>> for Geometry {
    fn from(l: &geo_types::Line<f64>) -> Self {
        Geometry::LineString(vec![coord(l.start), coord(l.end)])
    }
}

impl From<&geo_types::LineString<f64>> for Geometry {
    fn from(l: &geo_types::LineString<f64>) -> Self {
        Geometry::LineString(line_string(l))
    }
}

impl From<&geo_types::MultiLineString<f64>> for Geometry {
    fn from(l: &geo_types::MultiLineString<f64>) -> Self {
        Geometry::MultiLineString(l.0.iter().map(line_string).collect())
    }
}

impl From<&geo_types::Polygon<f64>> for Polygon {
    fn from(p: &geo_types::Polygon<f64>) -> Self {
        polygon(p)
    }
}

impl From<&geo_types::Polygon<f64>> for Geometry {
    fn from(p: &geo_types::Polygon<f64>) -> Self {
        Geometry::Polygon(polygon(p))
    }
}

impl From<&geo_types::MultiPolygon<f64>> for Geometry {
    fn from(p: &geo_types::MultiPolygon<f64>) -> Self {
        Geometry::MultiPolygon(p.0.iter().map(polygon).collect())
    }
}

impl From<&geo_types::Rect<f64>> for Geometry {
    fn from(r: &geo_types::Rect<f64>) -> Self {
        Geometry::Polygon(polygon(&r.to_polygon()))
    }
}

impl From<&geo_types::Triangle<f64>> for Geometry {
    fn from(t: &geo_types::Triangle<f64>) -> Self {
        Geometry::Polygon(polygon(&t.to_polygon()))
    }
}

impl TryFrom<&geo_types::Geometry<f64>> for Geometry {
    type Error = GeometryCollectionError;

    fn try_from(g: &geo_types::Geometry<f64>) -> Result<Self, Self::Error> {
        use geo_types::Geometry as G;
        Ok(match g {
            G::Point(p) => p.into(),
            G::MultiPoint(p) => p.into(),
            G::Line(l) => l.into(),
            G::LineString(l) => l.into(),
            G::MultiLineString(l) => l.into(),
            G::Polygon(p) => p.into(),
            G::MultiPolygon(p) => p.into(),
            G::Rect(r) => r.into(),
            G::Triangle(t) => t.into(),
            G::GeometryCollection(_) => return Err(GeometryCollectionError),
        })
    }
}

impl TryFrom<geo_types::Geometry<f64>> for Geometry {
    type Error = GeometryCollectionError;

    fn try_from(g: geo_types::Geometry<f64>) -> Result<Self, Self::Error> {
        Geometry::try_from(&g)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{burn, BurnOptions, GridSpec};

    #[test]
    fn polygon_round_trip_burns() {
        let ext = geo_types::LineString::from(vec![(2.5, 4.5), (6.5, 4.5), (6.5, 8.5), (2.5, 8.5)]);
        let p = geo_types::Polygon::new(ext, vec![]);
        let g: Geometry = (&p).into();
        assert_eq!(g.polygons()[0].rings[0].len(), 5, "geo-types closes rings");
        let grid = GridSpec::new(0., 0., 10., 10., 10, 10);
        let r = burn(&[g], &grid, BurnOptions::default()).unwrap();
        assert!((r.covered_cells() - 16.0).abs() < 1e-5);
    }

    #[test]
    fn collection_rejected() {
        let gc = geo_types::Geometry::GeometryCollection(geo_types::GeometryCollection::default());
        assert_eq!(Geometry::try_from(&gc), Err(GeometryCollectionError));
        let rect = geo_types::Geometry::Rect(geo_types::Rect::new((0., 0.), (1., 2.)));
        let g = Geometry::try_from(&rect).unwrap();
        assert_eq!(g.polygons().len(), 1);
    }
}
