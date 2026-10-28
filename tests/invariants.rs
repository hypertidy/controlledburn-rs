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

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon {
    Polygon::new(vec![[(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)].into_iter().map(Coord::from).collect()])
}

/// Coalesced runs are the parity runs passed through `coalesce_runs`, in
/// both modes, whatever the order of a multipolygon's parts.
fn assert_coalesced_matches_parity(g: &Geometry, grid: &GridSpec) {
    for mode in [BurnMode::Coverage, BurnMode::Approx] {
        let got = burn(std::slice::from_ref(g), grid, BurnOptions::default().with_mode(mode)).unwrap();
        let mut want = burn(std::slice::from_ref(g), grid, BurnOptions::parity(mode)).unwrap().runs;
        coalesce_runs(&mut want);
        assert_eq!(got.runs, want, "{mode:?}");
    }
}

#[test]
fn multipolygon_later_part_left_of_earlier_on_shared_row() {
    // Part 1 ends on row 5 at cols 6..9; part 2 starts on row 5 at cols
    // 1..4. Merging part 2's run into part 1's used to drop it.
    let g = Geometry::MultiPolygon(vec![rect(6., 4., 9., 7.), rect(1., 1., 4., 5.)]);
    let r = burn(std::slice::from_ref(&g), &G10, BurnOptions::default()).unwrap();
    assert_eq!(area(&r, &G10), 21.0);
    assert_coalesced_matches_parity(&g, &G10);
}

#[test]
fn multipolygon_touching_parts_coalesce_across_parts() {
    // The parts touch along x = 5 on row 5: one run, cols 2..8.
    let g = Geometry::MultiPolygon(vec![rect(5., 4., 8., 7.), rect(2., 2., 5., 5.)]);
    let r = burn(std::slice::from_ref(&g), &G10, BurnOptions::default()).unwrap();
    assert_eq!(area(&r, &G10), 18.0);
    let row5: Vec<_> = r.runs.iter().filter(|x| x.row == 5).copied().collect();
    assert_eq!(row5, vec![GridRun { row: 5, col_start: 2, col_end: 8, id: 0 }]);
    assert_coalesced_matches_parity(&g, &G10);
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
fn invalid_grids_are_rejected() {
    // Each of these used to pass validate(): inf and NaN bounds are caught
    // by the finiteness check; the last has finite bounds whose span
    // overflows f64, which gave an infinite cell size and an empty burn.
    let g = square(1., 1., 3., 3.);
    let bad = [
        GridSpec::new(0., 0., f64::INFINITY, 4., 4, 4),
        GridSpec::new(f64::NEG_INFINITY, 0., 4., 4., 4, 4),
        GridSpec::new(0., 0., f64::NAN, 4., 4, 4),
        GridSpec::new(0., f64::NAN, 4., 4., 4, 4),
        GridSpec::new(-1e308, 0., 1e308, 4., 4, 4),
        GridSpec::new(0., -1e308, 4., 1e308, 4, 4),
    ];
    for gs in bad {
        assert!(gs.validate().is_err(), "{gs:?} should be invalid");
        assert!(burn(std::slice::from_ref(&g), &gs, BurnOptions::default()).is_err(), "{gs:?}");
    }
    // a huge but representable span is still fine
    assert!(GridSpec::new(-1e307, 0., 1e307, 4., 4, 4).validate().is_ok());
}

/// Dense cover from the tables and from `CoverageWindow` (exact), for a
/// single ring; returns the worst cell as (row, col, burned, exact).
fn worst_cell_vs_window(ring: &[(f64, f64)], g: &GridSpec) -> (usize, usize, f64, f64) {
    let e = g.extent();
    let (nc, nr) = (g.ncol as usize, g.nrow as usize);
    let (dx, dy) = (g.dx(), g.dy());
    let n = ring.len() - 1;
    let xs: Vec<f64> = ring[..n].iter().map(|p| (p.0 - e.xmin) / dx).collect();
    let ys: Vec<f64> = ring[..n].iter().map(|p| (e.ymax - p.1) / dy).collect();
    let mut w = CoverageWindow::default();
    w.reset(0, 0, nc, nr);
    w.add_ring(&xs, &ys, 1.0);
    let exact = w.resolve().to_vec();
    let r = burn(&[poly(vec![ring.to_vec()])], g, BurnOptions::default()).unwrap();
    assert!(r.notes.is_empty(), "{:?}", r.notes);
    let mut m = vec![0.0; nc * nr];
    for x in &r.runs {
        for c in x.col_start..x.col_end {
            m[x.row as usize * nc + c as usize] += 1.0;
        }
    }
    for x in &r.edges {
        m[x.row as usize * nc + x.col as usize] += x.fraction as f64;
    }
    let i = (0..nc * nr).max_by(|&a, &b| (m[a] - exact[a]).abs().total_cmp(&(m[b] - exact[b]).abs())).unwrap();
    (i / nc, i % nc, m[i], exact[i])
}

fn assert_exact_every_start(ring: &[(f64, f64)], g: &GridSpec) {
    let n = ring.len() - 1;
    for reverse in [false, true] {
        for k in 0..n {
            let mut r: Vec<(f64, f64)> = (0..n).map(|i| ring[(i + k) % n]).collect();
            if reverse {
                r.reverse();
            }
            r.push(r[0]);
            let (row, col, got, want) = worst_cell_vs_window(&r, g);
            assert!(
                (got - want).abs() <= 2e-6,
                "start {k}, reversed {reverse}: row {row} col {col} burned {got} vs exact {want}"
            );
        }
    }
}

#[test]
fn edge_exactly_through_a_cell_corner_is_not_overcounted() {
    // controlledburn-rs#6 / controlledburn#25, from the property suite.
    // Each ring has an edge through a cell corner. The crossing came out an
    // ulp or two off the corner, the next cell got a sliver traversal along
    // its wall where the ring only touches the corner, and the chain chaser
    // read the sliver as a slab of the cell (or all of it).

    // (12.6, 1.4) -> (13.65, 0.35) passes through the top-right corner
    // (13.3, 0.7) of row 12 col 18: burned 1.0, exact 0.0943.
    assert_exact_every_start(
        &[
            (13.49307281327689, 0.280951025273552),
            (13.774064548591133, 0.1953615710514231),
            (14.150447247739818, 0.19536157105142304),
            (14.454947247739819, 0.2945346261241548),
            (14.571255898165475, 0.455),
            (14.454947247739819, 0.6154653738758452),
            (17.207761779641103, 4.932666216410308),
            (12.6, 1.4000000000000004),
            (13.649999999999999, 0.34999999999999964),
            (6.0614012054372814, 0.45500000000000046),
            (13.49307281327689, 0.280951025273552),
        ],
        &GridSpec::new(0.0, 0.0, 20.299999999999997, 9.1, 29, 13),
    );
    // (8.3, 10.4) -> (12.05, 6.65) passes through the corner (10.8, 7.9)
    // of row 4 col 4: burned 0.739, exact 0.455.
    assert_exact_every_start(
        &[
            (2.9404541633227357, 6.235726447802992),
            (7.140490483291759, -0.6880522184656641),
            (12.79922599805937, 6.649999999999999),
            (13.222799970968294, 7.999704574445243),
            (12.79922599805937, 8.351528741555857),
            (8.3, 10.399999999999999),
            (12.05, 6.649999999999999),
            (2.9404541633227357, 6.235726447802992),
        ],
        &GridSpec::new(-1.7, 2.9, 35.8, 17.9, 15, 6),
    );
    // (0.371, 0.8) -> (0.229, 0.657) passes through the corner (0.3, 0.729)
    // of row 3 col 0 on the grid's left edge: burned 0.839 or 1.0
    // depending on the start vertex, exact 0.339.
    assert_exact_every_start(
        &[
            (0.22857142857142856, 0.657142857142857),
            (0.3, 0.7285714285714284),
            (0.5385327833222915, 0.5669218770762636),
            (0.3826046772369731, 0.8968780211084503),
            (0.44285714285714284, 0.9150607669103901),
            (0.6141351860119378, 1.2110801232422874),
            (0.3, 0.9407750526246758),
            (0.37142857142857144, 0.7999999999999999),
            (0.22857142857142856, 0.657142857142857),
        ],
        &GridSpec::new(0.3, 0.3, 1.0142857142857142, 1.157142857142857, 5, 6),
    );
}

#[test]
fn near_coincident_edges_of_a_sliver_cover_nothing() {
    // controlledburn-rs#6, second case: three almost collinear lattice
    // points, area 4e-13. At 2x resolution its two long edges cross cell
    // walls an ulp apart in the wrong order, and the chain chaser returned
    // the complement of the sliver: 7 whole cells.
    let ring = [
        (1008.7666666666667, -2000.7),
        (1003.7666666666667, -1999.7),
        (1002.1, -1999.3666666666666),
        (1008.7666666666667, -2000.7),
    ];
    for k in [1u32, 2, 3, 4] {
        let g = GridSpec::new(1000.1, -2000.7, 1009.7666666666667, -1998.0333333333333, 29 * k, 8 * k);
        let r = burn(&[poly(vec![ring.to_vec()])], &g, BurnOptions::default()).unwrap();
        assert!(r.covered_cells() < 1e-5, "{k}x: covered {} cells", r.covered_cells());
        assert_exact_every_start(&ring, &g);
    }
    // Approx mode is the centre rule, exactly: the centre of row 4 col 8 is
    // strictly inside this ring (checked in rational arithmetic) and the
    // centre of row 6 col 18 lies on its edge, inside by the left-inclusive
    // rule. Both cells burn; see the CHANGELOG for why that stays.
    let g = GridSpec::new(1000.1, -2000.7, 1009.7666666666667, -1998.0333333333333, 29, 8);
    let r = burn(&[poly(vec![ring.to_vec()])], &g, BurnOptions::approx()).unwrap();
    let cells: Vec<(i32, i32)> =
        r.runs.iter().flat_map(|x| (x.col_start..x.col_end).map(move |c| (x.row, c))).collect();
    assert_eq!(cells, vec![(4, 8), (6, 18)]);
}

#[test]
fn approx_centre_tie_is_decided_exactly() {
    // From the property suite: the ring runs out along the line y = x - 3001.1
    // and back along it, and both edges pass the centre of row 13 col 2,
    // (1000.85, -2000.25), within 1.2e-13 on its right (exact rational
    // arithmetic). One rounded intercept came out exactly 1000.85, counted
    // as left of the centre by the left-inclusive rule, and the cell burned
    // with nothing in it.
    let g = GridSpec::new(1000.1, -2000.7, 1002.8000000000001, -1996.2, 9, 15);
    let ring = vec![
        (1000.4, -2000.7),
        (1001.3000000000001, -1999.8000000000002),
        (1001.293897355882, -1999.95),
        (1001.3511730051581, -1999.8000000000002),
        (1002.0812267040778, -1999.6101141954598),
        (1002.398382455067, -1997.7693056968183),
        (1001.293897355882, -1999.4751141954598),
        (1001.45, -1999.65),
        (1000.4, -2000.7),
    ];
    for reverse in [false, true] {
        let mut r = ring.clone();
        if reverse {
            r.reverse();
        }
        let res = burn(&[poly(vec![r])], &g, BurnOptions::approx()).unwrap();
        assert!(!res.runs.iter().any(|x| x.row == 13 && x.col_start <= 2 && x.col_end > 2), "{:?}", res.runs);
    }
}

#[test]
fn spike_to_a_corner_and_walls_traced_both_ways_stay_exact() {
    // More of controlledburn-rs#6, from the property suite at 100,000 cases.

    // A zero-width spike runs up to the top-left corner (0.2, 0.35) of row
    // 25 col 4 and straight back, returning an ulp from where it entered.
    // Its coverage came out 3e-16 instead of 0; that record blocked the
    // interior fill and the cell, fully covered, burned 0.
    assert_exact_every_start(
        &[
            (0.09000464930083663, 0.2951258728043881),
            (0.225, 0.2749999999999999),
            (0.2, 0.3500000000000001),
            (0.275, 0.125),
            (0.25, 0.44999999999999996),
            (0.2279343842256219, 0.4548286815756889),
            (0.22617971754901423, 0.4856824868406428),
            (0.22173674555812387, 0.5020994537202749),
            (0.2166843842256219, 0.49639790095734193),
            (0.09000464930083663, 0.2951258728043881),
        ],
        &GridSpec::new(0.0, 0.0, 0.25, 1.6, 5, 32),
    );
    // A triangle hanging below the grid traces part of the grid's bottom
    // edge out and back. From this start vertex the walk put one direction
    // in the bottom row and the other in the padding row below it; each
    // bottom-row cell then had a single traversal along its wall, which the
    // single-traversal formula read as the whole cell. Row 4 cols 17 and 18
    // burned 1, exact 0.
    assert_exact_every_start(
        &[
            (2.837265806357478, 0.30000000000000004),
            (3.014285714285714, 0.30000000000000004),
            (3.0402714541578892, 0.22857142857142854),
            (2.9374143113007465, 0.3),
            (2.837265806357478, 0.30000000000000004),
        ],
        &GridSpec::new(0.3, 0.3, 3.728571428571428, 1.0142857142857142, 24, 5),
    );
    // From CI on #17: an edge along lattice diagonals through the grid's
    // top-right corner. Row 2 col 9 burned as fully covered (exact 0.146),
    // so approx_agrees_with_coverage reported it missing from approx.
    let ring = [
        (-178.8212054697768, -87.7088728160751),
        (-178.7142857142857, -87.71428571428572),
        (-178.28571428571428, -87.28571428571429),
        (-178.64285714285714, -87.64285714285715),
        (-178.8212054697768, -87.61314778927319),
        (-178.8212054697768, -87.7088728160751),
    ];
    let g = GridSpec::new(-180.0, -90.0, -178.42857142857142, -87.28571428571429, 11, 19);
    assert_exact_every_start(&ring, &g);
    // Its centre is the vertex (-178.643, -87.643): inside by the
    // left-inclusive rule, so approx burns it.
    let r = burn(&[poly(vec![ring.to_vec()])], &g, BurnOptions::approx()).unwrap();
    assert!(r.runs.iter().any(|x| x.row == 2 && x.col_start <= 9 && x.col_end > 9));
}
