//! Port of the C++ core's test_burn.cpp plus the behaviours that differ
//! from the C++ by decision (ring closure, non-finite skipping, run
//! coalescing).

mod common;

use common::*;
use controlledburn::*;

fn poly(rings: Vec<Vec<(f64, f64)>>) -> Geometry {
    Geometry::Polygon(Polygon::new(rings.into_iter().map(|r| r.into_iter().map(Coord::from).collect()).collect()))
}

fn near(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * b.abs().max(1.0)
}

const G10: GridSpec = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);

#[test]
fn aligned_rectangle() {
    let r = burn(&[square(2., 4., 6., 8.)], &G10, BurnOptions::default()).unwrap();
    assert!(r.edges.is_empty());
    assert_eq!(r.runs.len(), 4, "coalesced: one run per row");
    assert_eq!(area(&r, &G10), 16.0);
    let p = burn(&[square(2., 4., 6., 8.)], &G10, BurnOptions::parity(BurnMode::Coverage)).unwrap();
    assert_eq!(p.runs.len(), 14, "C++ emits per-record runs");
    assert_eq!(area(&p, &G10), 16.0);
}

#[test]
fn offset_rectangle_conserves_area() {
    let r = burn(&[square(2.5, 4.5, 6.5, 8.5)], &G10, BurnOptions::default()).unwrap();
    assert!(!r.edges.is_empty());
    assert!(r.edges.iter().all(|e| e.fraction > 0.0 && e.fraction < 1.0));
    assert!(near(area(&r, &G10), 16.0, 1e-6));
}

#[test]
fn triangle_on_awkward_grid() {
    let g = GridSpec::new(0., 0., 100., 100., 37, 41);
    let t = poly(vec![vec![(13.3, 17.7), (88.1, 22.4), (41.9, 79.2), (13.3, 17.7)]]);
    let r = burn(std::slice::from_ref(&t), &g, BurnOptions::default()).unwrap();
    let exact = signed_area(&t.polygons()[0].rings[0]).abs();
    assert!(near(area(&r, &g), exact, 1e-4), "{} vs {exact}", area(&r, &g));
}

#[test]
fn holes_either_orientation() {
    let g = GridSpec::new(0., 0., 10., 10., 20, 20);
    for hole in
        [vec![(3., 3.), (7., 3.), (7., 7.), (3., 7.), (3., 3.)], vec![(3., 3.), (3., 7.), (7., 7.), (7., 3.), (3., 3.)]]
    {
        let p = poly(vec![vec![(1., 1.), (9., 1.), (9., 9.), (1., 9.), (1., 1.)], hole]);
        let r = burn(&[p], &g, BurnOptions::default()).unwrap();
        assert!(near(area(&r, &g), 48.0, 1e-6));
    }
}

#[test]
fn multipolygon_components_independent() {
    let a = Polygon::new(vec![vec![(1., 1.), (3., 1.), (3., 3.), (1., 3.), (1., 1.)]
        .into_iter()
        .map(Coord::from)
        .collect()]);
    let b = Polygon::new(vec![vec![(5.5, 5.5), (8.5, 5.5), (8.5, 8.5), (5.5, 8.5), (5.5, 5.5)]
        .into_iter()
        .map(Coord::from)
        .collect()]);
    let r = burn(&[Geometry::MultiPolygon(vec![a, b])], &G10, BurnOptions::default()).unwrap();
    assert!(near(area(&r, &G10), 13.0, 1e-6));
    assert!(r.runs.iter().chain(std::iter::empty()).all(|x| x.id == 1));
}

#[test]
fn beyond_extent_and_straddle() {
    let r = burn(&[square(-5., -5., 15., 15.)], &G10, BurnOptions::default()).unwrap();
    assert!(r.edges.is_empty());
    assert_eq!(area(&r, &G10), 100.0);
    assert_eq!(r.runs.len(), 10);
    let r = burn(&[square(7.5, 7.5, 12.5, 12.5)], &G10, BurnOptions::default()).unwrap();
    assert!(near(area(&r, &G10), 6.25, 1e-6));
    let r = burn(&[square(20., 20., 30., 30.)], &G10, BurnOptions::default()).unwrap();
    assert!(r.is_empty());
}

#[test]
fn line_length_conserved() {
    let l = Geometry::LineString(vec![Coord::new(0.5, 0.5), Coord::new(9.5, 7.5)]);
    let r = burn(&[l], &G10, BurnOptions::default()).unwrap();
    let want = (81.0f64 + 49.0).sqrt();
    assert!(near(r.line_length(), want, 1e-5));
    assert!(r.lines.iter().all(|c| c.length > 0.0));
    // Lines are unaffected by mode.
    let l = Geometry::LineString(vec![Coord::new(0.5, 0.5), Coord::new(9.5, 7.5)]);
    let r2 = burn(&[l], &G10, BurnOptions::approx()).unwrap();
    assert_eq!(r.lines, r2.lines);
}

#[test]
fn points_bin_and_drop() {
    let p = Geometry::MultiPoint(vec![
        Coord::new(0.5, 9.5),
        Coord::new(10.0, 0.0),
        Coord::new(10.5, 5.0),
        Coord::new(f64::NAN, 1.0),
    ]);
    // NaN inside a MultiPoint makes the whole geometry non-finite: skipped with a note.
    let r = burn(&[p], &G10, BurnOptions::default()).unwrap();
    assert!(r.points.is_empty());
    assert_eq!(r.notes.len(), 1);
    let p = Geometry::MultiPoint(vec![Coord::new(0.5, 9.5), Coord::new(10.0, 0.0), Coord::new(10.5, 5.0)]);
    let r = burn(&[p], &G10, BurnOptions::default()).unwrap();
    assert_eq!(r.points, vec![GridPoint { row: 1, col: 1, id: 1 }, GridPoint { row: 10, col: 10, id: 1 }]);
}

#[test]
fn ids_follow_input_position_and_empty_skipped() {
    let geoms = vec![
        Geometry::MultiPolygon(vec![]),
        square(1., 1., 3., 3.),
        Geometry::LineString(vec![]),
        Geometry::Point(Coord::new(5.5, 5.5)),
    ];
    let r = burn(&geoms, &G10, BurnOptions::default()).unwrap();
    assert!(r.runs.iter().all(|x| x.id == 2));
    assert_eq!(r.points[0].id, 4);
    assert!(r.notes.is_empty());
}

#[test]
fn wkb_roundtrip_and_malformed() {
    // POLYGON ((2 4, 6 4, 6 8, 2 8, 2 4)) plus a truncated blob and a collection
    let good = hex("010300000001000000050000000000000000000040000000000000104000000000000018400000000000001040000000000000184000000000000020400000000000000040000000000000204000000000000000400000000000001040");
    let bad = &good[..20];
    let coll = [1u8, 7, 0, 0, 0, 0, 0, 0, 0];
    let r = burn_wkb([good.as_slice(), bad, &coll, &[]], &G10, BurnOptions::default()).unwrap();
    assert_eq!(area(&r, &G10), 16.0);
    assert_eq!(r.notes.len(), 2);
    assert_eq!(r.notes[0].geom_index, 2);
    assert!(r.notes[0].message.starts_with("failed to parse WKB: WKB truncated at byte"));
    assert_eq!(r.notes[1].geom_index, 3);
    assert!(r.notes[1].message.contains("GeometryCollection"));
}

#[test]
fn materialize_policies() {
    let r = burn(&[square(2.5, 4.5, 6.5, 8.5)], &G10, BurnOptions::default()).unwrap();
    let mut buf = vec![f64::NAN; 100];
    materialize(&r, &mut buf, 10, 10, None, &MaterializeOptions::default()).unwrap();
    let touched = buf.iter().filter(|v| !v.is_nan()).count();
    // Threshold 0.5: interior 3x3 plus the 12 side cells at 0.5; the 4
    // corner cells at 0.25 are excluded.
    assert_eq!(touched, 21);
    assert!(buf.iter().filter(|v| !v.is_nan()).all(|v| *v == 1.0));

    let mut buf = vec![f64::NAN; 100];
    let opts = MaterializeOptions { fn_: PixelFn::Sum, edge_policy: EdgePolicy::Fraction, threshold: 0.5 };
    materialize(&r, &mut buf, 10, 10, Some(&[2.0]), &opts).unwrap();
    let total: f64 = buf.iter().filter(|v| !v.is_nan()).sum();
    assert!(near(total, 32.0, 1e-6));

    assert!(materialize(&r, &mut buf, 10, 9, None, &MaterializeOptions::default()).is_err());
    assert_eq!(
        materialize(&r, &mut buf, 10, 10, Some(&[]), &MaterializeOptions::default()).unwrap_err(),
        BurnError::IdOutOfRange { id: 1, values: 0 }
    );
}

#[test]
fn degenerate_inputs() {
    let collinear = poly(vec![vec![(1., 1.), (5., 5.), (9., 9.), (1., 1.)]]);
    let two = poly(vec![vec![(1., 1.), (5., 5.)]]);
    let r = burn(&[collinear, two], &G10, BurnOptions::default()).unwrap();
    assert!(r.is_empty() && r.notes.is_empty());
    assert!(burn(&[], &GridSpec::new(0., 0., 0., 10., 1, 1), BurnOptions::default()).is_err());
}

#[test]
fn unclosed_ring_is_closed() {
    let open = poly(vec![vec![(2., 2.), (6., 2.), (6., 6.), (2., 6.)]]);
    for opts in [BurnOptions::coverage(), BurnOptions::approx()] {
        let r = burn(std::slice::from_ref(&open), &G10, opts).unwrap();
        assert_eq!(area(&r, &G10), 16.0);
        assert!(r.notes.is_empty());
    }
}

#[test]
fn non_finite_polygon_is_skipped_not_hung() {
    let nan = poly(vec![vec![(2., 2.), (f64::NAN, 2.), (5., 5.), (2., 5.), (2., 2.)]]);
    let inf = poly(vec![vec![(2., 2.), (f64::INFINITY, 2.), (5., 5.), (2., 5.), (2., 2.)]]);
    let r = burn(&[nan, inf, square(1., 1., 2., 2.)], &G10, BurnOptions::default()).unwrap();
    assert_eq!(r.notes.len(), 2);
    assert_eq!((r.notes[0].geom_index, r.notes[1].geom_index), (1, 2));
    assert!(r.runs.iter().all(|x| x.id == 3));
}

#[test]
fn approx_matches_centre_rule() {
    // Offset rectangle: cell centres at .5 fall on the boundary and are
    // inside by the left-inclusive rule on the left/bottom edge only.
    let r = burn(&[square(2.5, 4.5, 6.5, 8.5)], &G10, BurnOptions::approx()).unwrap();
    assert!(r.edges.is_empty());
    assert_eq!(r.run_cells(), 16);
    assert_eq!(r.runs.len(), 4);
    // A rectangle whose edges avoid all centres burns exactly the cells
    // whose centres it contains: x in {2.5..6.5}, y in {4.5..8.5} = 5 x 5.
    let r = burn(&[square(2.2, 4.2, 6.7, 8.7)], &G10, BurnOptions::approx()).unwrap();
    assert_eq!(r.run_cells(), 25);
}
