//! Input from any [`geo-traits`](https://docs.rs/geo-traits) geometry
//! (feature `geo-traits`): geo-types, the `wkb` crate's zero-copy reader,
//! GeoArrow arrays, and anything else that implements `GeometryTrait`.
//!
//! [`Geometry::from_geo_traits`] copies one geometry into this crate's
//! [`Geometry`]; [`burn_geo_traits`](crate::burn_geo_traits) burns a
//! sequence of them, as [`burn_wkb`](crate::burn_wkb) burns WKB.
//!
//! Only `x` and `y` are read (Z and M are skipped, as in WKB). An empty
//! point is a point with both coordinates NaN, WKB's `POINT EMPTY`; a
//! polygon with no exterior is empty. `Rect` and `Triangle` become
//! polygons with the rings `geo-types`' `to_polygon` gives, and `Line` a
//! two-point line string. A `GeometryCollection` is refused, as by the
//! `geo-types` conversion: mixed dimensions have no place in the
//! four-table output, so split it upstream.

use geo_traits::{
    CoordTrait, GeometryTrait, GeometryType, LineStringTrait, LineTrait, MultiLineStringTrait, MultiPointTrait,
    MultiPolygonTrait, PointTrait, PolygonTrait, RectTrait, TriangleTrait,
};

use crate::geometry::{Coord, CoordSeq, Geometry, GeometryCollectionError, Polygon};

#[inline]
fn coord<C: CoordTrait<T = f64>>(c: &C) -> Coord {
    Coord::new(c.x(), c.y())
}

fn point<P: PointTrait<T = f64>>(p: &P) -> Coord {
    p.coord().map_or(Coord::new(f64::NAN, f64::NAN), |c| coord(&c))
}

fn line_string<L: LineStringTrait<T = f64>>(l: &L) -> CoordSeq {
    l.coords().map(|c| coord(&c)).collect()
}

fn polygon<P: PolygonTrait<T = f64>>(p: &P) -> Polygon {
    let Some(exterior) = p.exterior() else { return Polygon { rings: Vec::new() } };
    let mut rings = Vec::with_capacity(1 + p.num_interiors());
    rings.push(line_string(&exterior));
    rings.extend(p.interiors().map(|r| line_string(&r)));
    Polygon { rings }
}

/// A closed ring through `cs`, as `geo-types` builds one.
fn ring(cs: &[Coord]) -> Polygon {
    let mut r = cs.to_vec();
    r.push(cs[0]);
    Polygon { rings: vec![r] }
}

impl Geometry {
    /// Copy a `geo-traits` geometry (feature `geo-traits`). See the
    /// [module](crate::traits) for the rules; the only failure is a
    /// `GeometryCollection`.
    ///
    /// ```
    /// use controlledburn::{burn, BurnOptions, Geometry, GridSpec};
    ///
    /// // anything that implements geo_traits::GeometryTrait<T = f64>
    /// let ring = vec![(2.5, 4.5), (6.5, 4.5), (6.5, 8.5), (2.5, 8.5)];
    /// let square = geo_types::Geometry::Polygon(geo_types::Polygon::new(ring.into(), vec![]));
    /// let g = Geometry::from_geo_traits(&square).unwrap();
    /// let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
    /// let r = burn(&[g], &grid, BurnOptions::default()).unwrap();
    /// assert!((r.covered_cells() - 16.0).abs() < 1e-5);
    /// ```
    pub fn from_geo_traits<G: GeometryTrait<T = f64>>(g: &G) -> Result<Geometry, GeometryCollectionError> {
        Ok(match g.as_type() {
            GeometryType::Point(p) => Geometry::Point(point(p)),
            GeometryType::MultiPoint(m) => Geometry::MultiPoint(m.points().map(|p| point(&p)).collect()),
            GeometryType::Line(l) => Geometry::LineString(vec![coord(&l.start()), coord(&l.end())]),
            GeometryType::LineString(l) => Geometry::LineString(line_string(l)),
            GeometryType::MultiLineString(m) => {
                Geometry::MultiLineString(m.line_strings().map(|l| line_string(&l)).collect())
            }
            GeometryType::Polygon(p) => Geometry::Polygon(polygon(p)),
            GeometryType::MultiPolygon(m) => Geometry::MultiPolygon(m.polygons().map(|p| polygon(&p)).collect()),
            GeometryType::Rect(r) => {
                let (lo, hi) = (coord(&r.min()), coord(&r.max()));
                Geometry::Polygon(ring(&[
                    Coord::new(hi.x, lo.y),
                    Coord::new(hi.x, hi.y),
                    Coord::new(lo.x, hi.y),
                    Coord::new(lo.x, lo.y),
                ]))
            }
            GeometryType::Triangle(t) => {
                Geometry::Polygon(ring(&[coord(&t.first()), coord(&t.second()), coord(&t.third())]))
            }
            GeometryType::GeometryCollection(_) => return Err(GeometryCollectionError),
        })
    }
}
