//! `burn_geo_traits` (feature `geo-traits`) against the other inputs: the
//! `wkb` crate's zero-copy reader over every WKB fixture gives what
//! `burn_wkb` gives on the same bytes, and geo-types geometries through
//! geo-traits give what the `geo-types` conversion gives, for random
//! geometries of every type.
#![cfg(feature = "geo-traits")]

mod common;

use common::hex;
use controlledburn::*;
use proptest::prelude::*;

fn tables(r: &BurnResult) -> (&[GridRun], &[GridEdge], &[GridLine], &[GridPoint]) {
    (&r.runs, &r.edges, &r.lines, &r.points)
}

fn csv_rows(text: &str) -> Vec<Vec<String>> {
    text.lines()
        .skip(1)
        .filter(|l| !l.is_empty())
        .map(|l| l.split(',').map(|f| f.trim_matches('"').to_string()).collect())
        .collect()
}

/// Every WKB fixture, with its grid: burnscars' 50 oracle cases and the 14
/// cross-language parity cases.
fn wkb_cases() -> Vec<(String, Vec<u8>, GridSpec)> {
    let mut out = Vec::new();
    let oracle = std::fs::read_to_string(common::fixtures_dir().join("oracle/cases.csv")).unwrap();
    for r in csv_rows(&oracle) {
        let f = |i: usize| r[i].parse::<f64>().unwrap();
        let grid = GridSpec::new(f(2), f(3), f(4), f(5), r[6].parse().unwrap(), r[7].parse().unwrap());
        out.push((r[0].clone(), hex(&r[1]), grid));
    }
    // the WKT column holds commas: take the fields from the right
    let parity = std::fs::read_to_string(common::fixtures_dir().join("geometries.csv")).unwrap();
    for r in csv_rows(&parity) {
        let n = r.len();
        let f = |i: usize| r[i].parse::<f64>().unwrap();
        let grid =
            GridSpec::new(f(n - 7), f(n - 6), f(n - 5), f(n - 4), r[n - 3].parse().unwrap(), r[n - 2].parse().unwrap());
        out.push((r[0].clone(), hex(&r[n - 1]), grid));
    }
    out
}

#[test]
fn wkb_reader_equals_burn_wkb() {
    let cases = wkb_cases();
    assert_eq!(cases.len(), 64);
    for (name, blob, grid) in &cases {
        for mode in [BurnMode::Coverage, BurnMode::Approx] {
            for merge in [true, false] {
                let opts = BurnOptions::default().with_mode(mode).with_merge_parts(merge).with_id_offset(3);
                let want = burn_wkb([blob.as_slice()], grid, opts).unwrap();
                let g = wkb::reader::read_wkb(blob).unwrap();
                let got = burn_geo_traits([&g], grid, opts).unwrap();
                assert_eq!(got, want, "{name} {mode:?} merge {merge}");
            }
        }
    }
}

#[test]
fn empty_z_and_collections() {
    let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
    let opts = BurnOptions::default();
    // hand-made ISO WKB: POINT EMPTY, POINT Z, MULTIPOINT (one empty),
    // POLYGON EMPTY, big-endian LINESTRING
    let le = |v: f64| v.to_le_bytes().to_vec();
    let point =
        |t: u32, cs: &[f64]| [vec![1], t.to_le_bytes().to_vec(), cs.iter().flat_map(|&c| le(c)).collect()].concat();
    let empty = point(1, &[f64::NAN, f64::NAN]);
    let pz = point(1001, &[2.5, 3.5, 99.0]);
    let mp = [vec![1], 4u32.to_le_bytes().to_vec(), 2u32.to_le_bytes().to_vec(), point(1, &[1.5, 1.5]), empty.clone()]
        .concat();
    let poly_empty = [vec![1], 3u32.to_le_bytes().to_vec(), 0u32.to_le_bytes().to_vec()].concat();
    let be_line = [
        vec![0],
        2u32.to_be_bytes().to_vec(),
        2u32.to_be_bytes().to_vec(),
        [0.5f64, 0.5, 9.5, 7.25].iter().flat_map(|c| c.to_be_bytes()).collect(),
    ]
    .concat();
    let blobs = [empty, pz, mp, poly_empty, be_line];
    let want = burn_wkb(blobs.iter().map(|b| b.as_slice()), &grid, opts).unwrap();
    let gs: Vec<_> = blobs.iter().map(|b| wkb::reader::read_wkb(b).unwrap()).collect();
    let got = burn_geo_traits(gs.iter(), &grid, opts).unwrap();
    assert_eq!(got, want);
    assert_eq!(got.points.len(), 2);
    assert!(got.notes.is_empty());

    // a collection is skipped with a note, the rest burned
    let gc =
        geo_types::Geometry::GeometryCollection(geo_types::GeometryCollection(vec![
            geo_types::Point::new(1.0, 1.0).into()
        ]));
    let p = geo_types::Geometry::Point(geo_types::Point::new(4.5, 4.5));
    let r = burn_geo_traits([&gc, &p], &grid, opts).unwrap();
    assert_eq!(r.points, vec![GridPoint { row: 5, col: 4, id: 1 }]);
    assert_eq!(r.notes.len(), 1);
    assert_eq!(r.notes[0].geom_index, 0);
    assert_eq!(Geometry::from_geo_traits(&gc), Err(GeometryCollectionError));
}

fn ring(cx: f64, cy: f64, r: f64, k: usize, wobble: &[f64]) -> geo_types::LineString<f64> {
    (0..k)
        .map(|i| {
            let a = i as f64 / k as f64 * std::f64::consts::TAU;
            let rr = r * (0.6 + 0.4 * wobble[i % wobble.len()]);
            (cx + rr * a.cos(), cy + rr * a.sin())
        })
        .collect::<Vec<_>>()
        .into()
}

proptest! {
    #[test]
    fn geo_types_through_traits_equals_the_conversion(
        cx in 0.0..17.0f64, cy in 0.0..11.0f64, r in 0.5..9.0f64,
        k in 3usize..12, wobble in prop::collection::vec(0.0..1.0f64, 12),
        approx in any::<bool>(), merge in any::<bool>(),
    ) {
        use geo_types::{Geometry as G, Line, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon, Rect, Triangle};
        let grid = GridSpec::new(0.0, 0.0, 17.0, 11.0, 13, 7);
        let shell = ring(cx, cy, r, k, &wobble);
        let mut hole = ring(cx, cy, r * 0.3, k, &wobble);
        hole.0.reverse();
        let geoms: Vec<G<f64>> = vec![
            G::Polygon(Polygon::new(shell.clone(), vec![hole])),
            G::MultiPolygon(MultiPolygon(vec![
                Polygon::new(ring(cy, cx * 0.6, r * 0.8, k, &wobble), vec![]),
                Polygon::new(ring(cx * 0.5, cy * 0.5, r * 0.4, k, &wobble), vec![]),
            ])),
            G::Rect(Rect::new((cx - r, cy - r * 0.5), (cx + r * 0.7, cy + r))),
            G::Triangle(Triangle::new((cx, cy).into(), (cx + r, cy + 0.3).into(), (cx - 0.2, cy + r).into())),
            G::LineString(shell.clone()),
            G::MultiLineString(MultiLineString(vec![shell, ring(cy, cx, r * 0.7, k, &wobble)])),
            G::Line(Line::new((cx, cy), (cy, cx))),
            G::MultiPoint(MultiPoint((0..k).map(|i| Point::new(cx + i as f64, cy - i as f64 * 0.5)).collect())),
            G::Point(Point::new(cx, cy)),
        ];
        let mode = if approx { BurnMode::Approx } else { BurnMode::Coverage };
        let opts = BurnOptions::default().with_mode(mode).with_merge_parts(merge);
        let converted: Vec<Geometry> = geoms.iter().map(|g| Geometry::try_from(g).unwrap()).collect();
        let want = burn(&converted, &grid, opts).unwrap();
        let got = burn_geo_traits(geoms.iter(), &grid, opts).unwrap();
        prop_assert_eq!(tables(&got), tables(&want));
        prop_assert_eq!(&got.notes, &want.notes);
        for (g, c) in geoms.iter().zip(&converted) {
            prop_assert_eq!(&Geometry::from_geo_traits(g).unwrap(), c);
        }
    }
}
