//! Minimal WKB reader for the geometry model.
//!
//! Supports ISO WKB and EWKB: both byte orders, 2-D coordinates with Z
//! and/or M ordinates skipped (ISO 1000/2000/3000 offsets and EWKB
//! 0x80000000 / 0x40000000 flags), EWKB embedded SRID skipped. Curved
//! types and GeometryCollection are rejected: linearise / split upstream.

use crate::error::WkbError;
use crate::geometry::{Coord, CoordSeq, Geometry, Polygon};

struct Cursor<'a> {
    data: &'a [u8],
    pos: usize,
    little: bool,
}

impl<'a> Cursor<'a> {
    fn new(data: &'a [u8]) -> Self {
        Cursor { data, pos: 0, little: true }
    }

    fn remaining(&self) -> usize {
        self.data.len() - self.pos
    }

    fn need(&self, n: usize) -> Result<(), WkbError> {
        if self.pos + n > self.data.len() {
            Err(WkbError::Truncated { at: self.pos })
        } else {
            Ok(())
        }
    }

    fn read_u8(&mut self) -> Result<u8, WkbError> {
        self.need(1)?;
        let v = self.data[self.pos];
        self.pos += 1;
        Ok(v)
    }

    fn read_u32(&mut self) -> Result<u32, WkbError> {
        self.need(4)?;
        let b: [u8; 4] = self.data[self.pos..self.pos + 4].try_into().unwrap();
        self.pos += 4;
        Ok(if self.little { u32::from_le_bytes(b) } else { u32::from_be_bytes(b) })
    }

    fn read_f64(&mut self) -> Result<f64, WkbError> {
        self.need(8)?;
        let b: [u8; 8] = self.data[self.pos..self.pos + 8].try_into().unwrap();
        self.pos += 8;
        Ok(if self.little { f64::from_le_bytes(b) } else { f64::from_be_bytes(b) })
    }

    fn skip_f64(&mut self, n: usize) -> Result<(), WkbError> {
        self.need(8 * n)?;
        self.pos += 8 * n;
        Ok(())
    }

    fn read_byte_order(&mut self) -> Result<(), WkbError> {
        let order = self.read_u8()?;
        match order {
            0 => self.little = false,
            1 => self.little = true,
            o => return Err(WkbError::BadByteOrder(o)),
        }
        Ok(())
    }
}

const EWKB_Z: u32 = 0x8000_0000;
const EWKB_M: u32 = 0x4000_0000;
const EWKB_SRID: u32 = 0x2000_0000;

struct TypeInfo {
    base: u32,
    extra_dims: usize,
    has_srid: bool,
}

fn decode_type(raw: u32) -> TypeInfo {
    let mut extra = 0;
    if raw & EWKB_Z != 0 {
        extra += 1;
    }
    if raw & EWKB_M != 0 {
        extra += 1;
    }
    let mut base = raw & 0x0FFF_FFFF;
    if (3000..4000).contains(&base) {
        extra += 2;
        base -= 3000;
    } else if (2000..3000).contains(&base) {
        extra += 1;
        base -= 2000;
    } else if (1000..2000).contains(&base) {
        extra += 1;
        base -= 1000;
    }
    TypeInfo { base, extra_dims: extra, has_srid: raw & EWKB_SRID != 0 }
}

fn read_header(cur: &mut Cursor<'_>) -> Result<TypeInfo, WkbError> {
    cur.read_byte_order()?;
    let t = decode_type(cur.read_u32()?);
    if t.has_srid {
        cur.read_u32()?;
    }
    Ok(t)
}

fn read_nested_header(cur: &mut Cursor<'_>, expected_base: u32) -> Result<TypeInfo, WkbError> {
    let t = read_header(cur)?;
    if t.base != expected_base {
        return Err(WkbError::MismatchedPart { found: t.base });
    }
    Ok(t)
}

fn read_point(cur: &mut Cursor<'_>, extra: usize) -> Result<Coord, WkbError> {
    let x = cur.read_f64()?;
    let y = cur.read_f64()?;
    if extra > 0 {
        cur.skip_f64(extra)?;
    }
    Ok(Coord { x, y })
}

fn read_coordseq(cur: &mut Cursor<'_>, extra: usize) -> Result<CoordSeq, WkbError> {
    let n = cur.read_u32()? as usize;
    let mut seq = Vec::with_capacity(n.min(cur.remaining() / 16));
    for _ in 0..n {
        seq.push(read_point(cur, extra)?);
    }
    Ok(seq)
}

fn read_polygon_body(cur: &mut Cursor<'_>, extra: usize) -> Result<Polygon, WkbError> {
    let n_rings = cur.read_u32()? as usize;
    let mut rings = Vec::with_capacity(n_rings);
    for _ in 0..n_rings {
        rings.push(read_coordseq(cur, extra)?);
    }
    Ok(Polygon { rings })
}

/// Parse a single WKB geometry.
pub fn parse_wkb(data: &[u8]) -> Result<Geometry, WkbError> {
    if data.len() < 5 {
        return Err(WkbError::TooShort);
    }
    let mut cur = Cursor::new(data);
    let t = read_header(&mut cur)?;

    match t.base {
        1 => Ok(Geometry::Point(read_point(&mut cur, t.extra_dims)?)),
        2 => Ok(Geometry::LineString(read_coordseq(&mut cur, t.extra_dims)?)),
        3 => Ok(Geometry::Polygon(read_polygon_body(&mut cur, t.extra_dims)?)),
        4 => {
            let n = cur.read_u32()? as usize;
            let mut pts = Vec::with_capacity(n);
            for _ in 0..n {
                let pt = read_nested_header(&mut cur, 1)?;
                pts.push(read_point(&mut cur, pt.extra_dims)?);
            }
            Ok(Geometry::MultiPoint(pts))
        }
        5 => {
            let n = cur.read_u32()? as usize;
            let mut lines = Vec::with_capacity(n);
            for _ in 0..n {
                let lt = read_nested_header(&mut cur, 2)?;
                lines.push(read_coordseq(&mut cur, lt.extra_dims)?);
            }
            Ok(Geometry::MultiLineString(lines))
        }
        6 => {
            let n = cur.read_u32()? as usize;
            let mut polys = Vec::with_capacity(n);
            for _ in 0..n {
                let pt = read_nested_header(&mut cur, 3)?;
                polys.push(read_polygon_body(&mut cur, pt.extra_dims)?);
            }
            Ok(Geometry::MultiPolygon(polys))
        }
        7 => Err(WkbError::GeometryCollection),
        other => Err(WkbError::Unsupported { type_code: other }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hex(s: &str) -> Vec<u8> {
        (0..s.len()).step_by(2).map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap()).collect()
    }

    #[test]
    fn polygon_le() {
        // POLYGON ((2 4, 6 4, 6 8, 2 8, 2 4))
        let b = hex("010300000001000000050000000000000000000040000000000000104000000000000018400000000000001040000000000000184000000000000020400000000000000040000000000000204000000000000000400000000000001040");
        let g = parse_wkb(&b).unwrap();
        match g {
            Geometry::Polygon(p) => {
                assert_eq!(p.rings.len(), 1);
                assert_eq!(p.rings[0][1], Coord::new(6., 4.));
            }
            _ => panic!(),
        }
    }

    #[test]
    fn point_be_and_ewkb_zm_srid() {
        // XDR POINT (1 2)
        let b = hex("00000000013ff00000000000004000000000000000");
        assert_eq!(parse_wkb(&b).unwrap(), Geometry::Point(Coord::new(1., 2.)));
        // EWKB POINT ZM with SRID (little endian): type = 1 | Z | M | SRID
        let mut v = vec![1u8];
        v.extend_from_slice(&(1u32 | EWKB_Z | EWKB_M | EWKB_SRID).to_le_bytes());
        v.extend_from_slice(&4326u32.to_le_bytes());
        for f in [3.0f64, 4.0, 5.0, 6.0] {
            v.extend_from_slice(&f.to_le_bytes());
        }
        assert_eq!(parse_wkb(&v).unwrap(), Geometry::Point(Coord::new(3., 4.)));
        // ISO POINT Z (1001)
        let mut v = vec![1u8];
        v.extend_from_slice(&1001u32.to_le_bytes());
        for f in [3.0f64, 4.0, 5.0] {
            v.extend_from_slice(&f.to_le_bytes());
        }
        assert_eq!(parse_wkb(&v).unwrap(), Geometry::Point(Coord::new(3., 4.)));
    }

    #[test]
    fn errors() {
        assert_eq!(parse_wkb(&[1, 1, 0]).unwrap_err(), WkbError::TooShort);
        assert_eq!(parse_wkb(&[2, 1, 0, 0, 0, 0]).unwrap_err(), WkbError::BadByteOrder(2));
        assert_eq!(parse_wkb(&[1, 7, 0, 0, 0, 0]).unwrap_err(), WkbError::GeometryCollection);
        assert_eq!(parse_wkb(&[1, 8, 0, 0, 0, 0]).unwrap_err(), WkbError::Unsupported { type_code: 8 });
        assert!(matches!(parse_wkb(&[1, 1, 0, 0, 0, 0, 0]).unwrap_err(), WkbError::Truncated { .. }));
        // MultiPoint with a LineString part
        let mut v = vec![1u8];
        v.extend_from_slice(&4u32.to_le_bytes());
        v.extend_from_slice(&1u32.to_le_bytes());
        v.push(1);
        v.extend_from_slice(&2u32.to_le_bytes());
        assert_eq!(parse_wkb(&v).unwrap_err(), WkbError::MismatchedPart { found: 2 });
    }
}
