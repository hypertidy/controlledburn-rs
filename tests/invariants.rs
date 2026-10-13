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
    assert!(r.runs.iter().all(|x| x.id == 0));
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
    assert_eq!(r.points, vec![GridPoint { row: 0, col: 0, id: 0 }, GridPoint { row: 9, col: 9, id: 0 }]);
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
    assert!(r.runs.iter().all(|x| x.id == 1));
    assert_eq!(r.points[0].id, 3);
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
    assert_eq!(r.notes[0].geom_index, 1);
    assert!(r.notes[0].message.starts_with("failed to parse WKB: WKB truncated at byte"));
    assert_eq!(r.notes[1].geom_index, 2);
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
    assert!(buf.iter().filter(|v| !v.is_nan()).all(|v| *v == 0.0), "burned value is the id, 0");

    let mut buf = vec![f64::NAN; 100];
    let opts = MaterializeOptions { fn_: PixelFn::Sum, edge_policy: EdgePolicy::Fraction, threshold: 0.5 };
    materialize(&r, &mut buf, 10, 10, Some(&[2.0]), &opts).unwrap();
    let total: f64 = buf.iter().filter(|v| !v.is_nan()).sum();
    assert!(near(total, 32.0, 1e-6));

    assert!(materialize(&r, &mut buf, 10, 9, None, &MaterializeOptions::default()).is_err());
    assert_eq!(
        materialize(&r, &mut buf, 10, 10, Some(&[]), &MaterializeOptions::default()).unwrap_err(),
        BurnError::IdOutOfRange { id: 0, values: 0 }
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
    assert_eq!((r.notes[0].geom_index, r.notes[1].geom_index), (0, 1));
    assert!(r.runs.iter().all(|x| x.id == 2));
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

#[test]
fn vertex_a_few_ulp_off_a_cell_edge_does_not_derail_the_walk() {
    // 0.3 parsed from input vs 3 * 0.1 computed from the grid differ by a
    // few ulp. Before lattice snapping the walk left the grid along the
    // extrapolated edge and ran until memory was exhausted (the same input
    // crashes exactextractr). Checked against the analytic triangle area.
    let g = GridSpec::new(0.0, 0.0, 10.0, 10.0, 100, 100);
    let tri = |a: (f64, f64), b: (f64, f64), c: (f64, f64)| {
        ((b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1)).abs() / 2.0
    };
    let cases = [
        (vec![(0.3, 0.3), (4.7, 0.3), (2.5, 6.2), (0.3, 0.3)], tri((0.3, 0.3), (4.7, 0.3), (2.5, 6.2))),
        (vec![(4.7, 0.3), (2.5, 6.2), (0.3, 0.3), (4.7, 0.3)], tri((0.3, 0.3), (4.7, 0.3), (2.5, 6.2))),
        (vec![(0.3, 0.3), (4.7, 0.3), (2.5, 6.25), (0.3, 0.3)], tri((0.3, 0.3), (4.7, 0.3), (2.5, 6.25))),
    ];
    for (ring, expected) in cases {
        let r = burn(&[poly(vec![ring])], &g, BurnOptions::default()).unwrap();
        assert!(!r.runs.is_empty());
        assert!(near(area(&r, &g), expected, 1e-6), "area {} vs {}", area(&r, &g), expected);
    }
}

#[test]
fn row_centre_crossing_on_a_column_wall_keeps_its_winding() {
    // Falkland Islands from Natural Earth 110m, shifted so every vertex
    // sits on a 0.05 multiple. At 20 columns (res 0.5) the first edge
    // crosses the row centre y = 1.25 at x = 1.0, exactly on a column wall;
    // with a strict crossing test neither adjacent cell registered it and
    // the row lost two interior cells. Exact at every resolution now.
    let ring = vec![
        (0.8, 1.15),
        (2.0, 1.75),
        (2.85, 1.5),
        (3.45, 1.9),
        (4.25, 1.45),
        (3.95, 1.1),
        (2.6, 0.8),
        (2.15, 1.15),
        (1.3, 0.7),
        (0.8, 1.15),
    ];
    let expected = 2.12875;
    for n in [20u32, 100, 200, 1000] {
        let g = GridSpec::new(0.0, 0.0, 10.0, 10.0, n, n);
        let r = burn(&[poly(vec![ring.clone()])], &g, BurnOptions::default()).unwrap();
        assert!(near(area(&r, &g), expected, 1e-6), "ncol {n}: area {} vs {}", area(&r, &g), expected);
    }
}

#[test]
fn approx_vertex_on_row_centre_does_not_fill_to_grid_edge() {
    // Found by the property suite (tests/properties.rs), shrunk by proptest.
    // The middle vertex sits on the centre line of row 24 as the generator
    // computed it (ymax - 24.5 / 3), one ulp below the same line recomputed
    // from the extent-derived dy. The row-range test admitted the downward
    // edge starting there; the sign test (y0 >= y_mid) called it upward;
    // the two intercepts then had the same sign and the row filled from
    // the polygon to the right edge of the grid.
    let g = GridSpec::new(-180.0, -90.0, -175.0, -80.66666666666667, 15, 28);
    let ring = vec![
        (-179.89654756364214, -88.83333333333334),
        (-180.94076559960791, -92.3680897750401),
        (-179.6, -88.518152273978),
        (-179.89654756364214, -88.83333333333334),
    ];
    let r = burn(&[poly(vec![ring])], &g, BurnOptions::approx()).unwrap();
    for run in &r.runs {
        assert!(run.col_end <= 2, "approx run spills past the polygon: {run:?}");
    }
    assert_eq!(r.runs.len(), 2);
}

#[test]
fn hole_boundary_cells_keep_the_part_outside_the_hole() {
    // A cell traversed only by a hole ring carried a negative hole fraction
    // and no exterior fraction, so the sweep emitted nothing for it: every
    // polygon with a hole lost the ring of cells along the hole boundary.
    let outer = vec![(0.5, 0.5), (9.5, 0.5), (9.5, 9.5), (0.5, 9.5), (0.5, 0.5)];
    let hole = vec![(3.3, 3.3), (6.7, 3.3), (6.7, 6.7), (3.3, 6.7), (3.3, 3.3)];
    let truth = 81.0 - 3.4 * 3.4;
    for n in [10u32, 20, 50, 100] {
        let g = GridSpec::new(0.0, 0.0, 10.0, 10.0, n, n);
        let r = burn(&[poly(vec![outer.clone(), hole.clone()])], &g, BurnOptions::default()).unwrap();
        assert!(near(area(&r, &g), truth, 1e-6), "ncol {n}: area {} vs {truth}", area(&r, &g));
    }
    // at 10 columns, cell (row 6, col 3) holds the hole's corner: covered
    // except 0.7 x 0.7 of it
    let g = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
    let r = burn(&[poly(vec![outer, hole])], &g, BurnOptions::default()).unwrap();
    let e: Vec<_> = r.edges.iter().filter(|e| e.row == 6 && e.col == 3).collect();
    assert_eq!(e.len(), 1);
    assert!((e[0].fraction as f64 - 0.51).abs() < 1e-6, "fraction {}", e[0].fraction);
}
