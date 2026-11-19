//! The order of the records, which `BurnResult` documents as part of the
//! contract: each table grouped by id in input order; within an id, runs,
//! edges and lines in row-major order (row, then column) with one record
//! per cell, and points in the order the geometry lists them. With
//! `merge_parts` off (and so in parity mode) that holds per part, and a
//! single-part geometry is still in row-major order. `merge_cells` leaves
//! a default burn as it is. `burn_par` and
//! `burn_wkb_par` equal `burn` record for record (tests in src/burn.rs),
//! and `burn_window`, `crop`, `slice_tiles` and `slice_tiles_par` keep the
//! order.

use controlledburn::{
    burn, burn_par, burn_window, burn_wkb, burn_wkb_par, crop, merge_cells, slice_tiles, slice_tiles_par, BurnMode,
    BurnOptions, BurnResult, Coord, Geometry, GridSpec, Polygon,
};
use proptest::prelude::*;

/// Ordered as documented, and no cell twice for an id.
fn check(r: &BurnResult) -> Result<(), String> {
    let strictly = |name: &str, keys: Vec<(i32, i32, i32)>| match keys.windows(2).find(|w| w[0] >= w[1]) {
        Some(w) => Err(format!("{name}: {:?} then {:?}", w[0], w[1])),
        None => Ok(()),
    };
    strictly("runs", r.runs.iter().map(|x| (x.id, x.row, x.col_start)).collect())?;
    strictly("edges", r.edges.iter().map(|x| (x.id, x.row, x.col)).collect())?;
    strictly("lines", r.lines.iter().map(|x| (x.id, x.row, x.col)).collect())?;
    if let Some(w) = r.points.windows(2).find(|w| w[0].id > w[1].id) {
        return Err(format!("points: id {} then {}", w[0].id, w[1].id));
    }
    // Runs of an id don't overlap, and no edge is in a run.
    if let Some(w) =
        r.runs.windows(2).find(|w| (w[0].id, w[0].row) == (w[1].id, w[1].row) && w[1].col_start < w[0].col_end)
    {
        return Err(format!("runs overlap: {:?} {:?}", w[0], w[1]));
    }
    for e in &r.edges {
        let i = r.runs.partition_point(|x| (x.id, x.row, x.col_start) <= (e.id, e.row, e.col));
        if i > 0 && (r.runs[i - 1].id, r.runs[i - 1].row) == (e.id, e.row) && e.col < r.runs[i - 1].col_end {
            return Err(format!("edge {e:?} in run {:?}", r.runs[i - 1]));
        }
    }
    Ok(())
}

fn coord() -> impl Strategy<Value = Coord> {
    (-2.0..32.0f64, -2.0..22.0f64).prop_map(|(x, y)| Coord::new(x, y))
}

/// A star-shaped ring around a centre: simple, either orientation.
fn star(cx: f64, cy: f64, radii: &[f64], ccw: bool) -> Vec<Coord> {
    let n = radii.len();
    let mut v: Vec<Coord> = radii
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let a = i as f64 / n as f64 * std::f64::consts::TAU;
            Coord::new(cx + r * a.cos(), cy + r * a.sin())
        })
        .collect();
    if !ccw {
        v.reverse();
    }
    v.push(v[0]);
    v
}

fn polygon() -> impl Strategy<Value = Polygon> {
    (coord(), prop::collection::vec(1.5..8.0f64, 3..14), any::<bool>(), prop::option::of(0.2..0.6f64)).prop_map(
        |(c, radii, ccw, hole)| {
            let mut rings = vec![star(c.x, c.y, &radii, ccw)];
            if let Some(k) = hole {
                // inside the outer ring: the shortest radius scaled down
                let r = radii.iter().cloned().fold(f64::INFINITY, f64::min) * k;
                rings.push(star(c.x, c.y, &[r; 6], !ccw));
            }
            Polygon::new(rings)
        },
    )
}

fn line() -> impl Strategy<Value = Vec<Coord>> {
    prop::collection::vec(coord(), 2..9)
}

fn geometry() -> impl Strategy<Value = Geometry> {
    prop_oneof![
        polygon().prop_map(Geometry::Polygon),
        prop::collection::vec(polygon(), 2..4).prop_map(Geometry::MultiPolygon),
        line().prop_map(Geometry::LineString),
        prop::collection::vec(line(), 2..4).prop_map(Geometry::MultiLineString),
        coord().prop_map(Geometry::Point),
        prop::collection::vec(coord(), 2..8).prop_map(Geometry::MultiPoint),
    ]
}

fn single_part(g: &Geometry) -> bool {
    !matches!(g, Geometry::MultiPolygon(_) | Geometry::MultiLineString(_) | Geometry::MultiPoint(_))
}

fn grid() -> impl Strategy<Value = GridSpec> {
    (8u32..40, 6u32..30).prop_map(|(nc, nr)| GridSpec::new(0.0, 0.0, 30.0, 20.0, nc, nr))
}

fn wkb(g: &Geometry) -> Vec<u8> {
    fn seq(o: &mut Vec<u8>, s: &[Coord]) {
        o.extend((s.len() as u32).to_le_bytes());
        s.iter().for_each(|c| o.extend(c.x.to_le_bytes().into_iter().chain(c.y.to_le_bytes())));
    }
    let mut o = vec![1u8];
    let (t, parts): (u32, Vec<Geometry>) = match g {
        Geometry::Point(c) => {
            o.extend(1u32.to_le_bytes());
            o.extend(c.x.to_le_bytes().into_iter().chain(c.y.to_le_bytes()));
            return o;
        }
        Geometry::LineString(s) => {
            o.extend(2u32.to_le_bytes());
            seq(&mut o, s);
            return o;
        }
        Geometry::Polygon(p) => {
            o.extend(3u32.to_le_bytes());
            o.extend((p.rings.len() as u32).to_le_bytes());
            p.rings.iter().for_each(|r| seq(&mut o, r));
            return o;
        }
        Geometry::MultiPoint(ps) => (4, ps.iter().map(|&c| Geometry::Point(c)).collect()),
        Geometry::MultiLineString(ls) => (5, ls.iter().map(|l| Geometry::LineString(l.clone())).collect()),
        Geometry::MultiPolygon(ps) => (6, ps.iter().map(|p| Geometry::Polygon(p.clone())).collect()),
    };
    o.extend(t.to_le_bytes());
    o.extend((parts.len() as u32).to_le_bytes());
    parts.iter().for_each(|p| o.extend(wkb(p)));
    o
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 400, ..ProptestConfig::default() })]

    #[test]
    fn records_are_in_contract_order(
        geoms in prop::collection::vec(geometry(), 1..7),
        grid in grid(),
        approx in any::<bool>(),
        coalesce in any::<bool>(),
        (r0, nr, c0, nc) in (0u32..6, 1u32..12, 0u32..8, 1u32..16),
        (tx, ty) in (3u32..12, 2u32..9),
        target in (0.0..15.0f64, 0.0..10.0f64, 1.0..15.0f64, 1.0..10.0f64),
    ) {
        let mode = if approx { BurnMode::Approx } else { BurnMode::Coverage };
        let opts = BurnOptions::default().with_mode(mode).with_coalesce_runs(coalesce);
        let blobs: Vec<Vec<u8>> = geoms.iter().map(wkb).collect();

        let r = burn(&geoms, &grid, opts).unwrap();
        check(&r).map_err(|e| TestCaseError::fail(format!("burn: {e}")))?;
        let w = burn_wkb(blobs.iter().map(|b| b.as_slice()), &grid, opts).unwrap();
        prop_assert_eq!(&w, &r);
        prop_assert_eq!(&burn_par(&geoms, &grid, opts, 3).unwrap(), &r);
        prop_assert_eq!(&burn_wkb_par(&blobs, &grid, opts, 3).unwrap(), &r);
        if coalesce {
            // already merged and in merge_cells' order
            let mut m = r.clone();
            merge_cells(&mut m);
            prop_assert_eq!(&m, &r);
        }

        let (rows, cols) = (r0..(r0 + nr).min(grid.nrow), c0..(c0 + nc).min(grid.ncol));
        let win = burn_window(&geoms, &grid, rows, cols, opts).unwrap();
        check(&win).map_err(|e| TestCaseError::fail(format!("burn_window: {e}")))?;

        let (x, y, w, h) = target;
        if let Some((c, _)) = crop(&r, &grid, [x, y, x + w, y + h]) {
            check(&c).map_err(|e| TestCaseError::fail(format!("crop: {e}")))?;
        }
        let tiles = slice_tiles(&r, &grid, tx, ty);
        for t in &tiles {
            check(&t.result).map_err(|e| TestCaseError::fail(format!("slice_tiles: {e}")))?;
        }
        let par = slice_tiles_par(&r, &grid, tx, ty, 3);
        prop_assert_eq!(par.len(), tiles.len());
        for (a, b) in par.iter().zip(&tiles) {
            prop_assert_eq!(&a.result, &b.result);
        }
    }

    /// With merging off, as in parity mode, single-part geometries are
    /// still in row-major order; multi-part ones are per part.
    #[test]
    fn parity_single_part_order(
        geoms in prop::collection::vec(geometry().prop_filter("single part", single_part), 1..7),
        grid in grid(),
        approx in any::<bool>(),
    ) {
        let mode = if approx { BurnMode::Approx } else { BurnMode::Coverage };
        let r = burn(&geoms, &grid, BurnOptions::parity(mode)).unwrap();
        // Parity runs are not coalesced, so adjacent runs may touch; they
        // still don't overlap and are in order.
        check(&r).map_err(|e| TestCaseError::fail(format!("parity: {e}")))?;
    }
}

/// The figure in the book's contract chapter: two polygons, the order of
/// their records, and a cell they share.
#[test]
fn two_polygons_of_the_figure() {
    let ring = |v: &[(f64, f64)]| v.iter().map(|&(x, y)| Coord::new(x, y)).collect::<Vec<_>>();
    let geoms = [
        Geometry::Polygon(Polygon::new(vec![
            ring(&[(0.6, 0.7), (5.3, 0.7), (7.4, 5.3), (0.6, 5.3), (0.6, 0.7)]),
            ring(&[(2.4, 2.4), (2.4, 3.6), (4.6, 3.6), (4.6, 2.4), (2.4, 2.4)]),
        ])),
        Geometry::Polygon(Polygon::new(vec![ring(&[(5.3, 0.7), (9.4, 0.7), (9.4, 5.3), (7.4, 5.3), (5.3, 0.7)])])),
    ];
    let grid = GridSpec::new(0.0, 0.0, 10.0, 6.0, 10, 6);
    let r = burn(&geoms, &grid, BurnOptions::default()).unwrap();
    check(&r).unwrap();
    assert_eq!((r.runs.len(), r.edges.len()), (9, 48));
    let shared: Vec<_> = r.edges.iter().filter(|e| (e.row, e.col) == (2, 6)).collect();
    assert_eq!(shared.iter().map(|e| e.id).collect::<Vec<_>>(), vec![0, 1]);
    assert!((shared[0].fraction + shared[1].fraction - 1.0).abs() < 1e-6);
}
