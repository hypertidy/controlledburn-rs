//! One record per (row, col, id) for multi-part geometries, notes for the
//! rings and lines the engines skip, and `id_offset`.

use std::collections::HashMap;

use controlledburn::{
    burn, burn_wkb, BurnError, BurnMode, BurnOptions, BurnResult, Coord, Geometry, GridSpec, Polygon,
};

fn g10() -> GridSpec {
    GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10)
}

fn ring(v: &[(f64, f64)]) -> Vec<Coord> {
    v.iter().map(|&(x, y)| Coord::new(x, y)).collect()
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon {
    Polygon::new(vec![ring(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)])])
}

fn multi(parts: Vec<Polygon>) -> Geometry {
    Geometry::MultiPolygon(parts)
}

/// Records per cell and the cell's value (run 1, edge its fraction),
/// keyed by (row, col, id).
fn cells(r: &BurnResult) -> HashMap<(i32, i32, i32), (usize, f64)> {
    let mut m: HashMap<(i32, i32, i32), (usize, f64)> = HashMap::new();
    for x in &r.runs {
        for c in x.col_start..x.col_end {
            let e = m.entry((x.row, c, x.id)).or_default();
            e.0 += 1;
            e.1 += 1.0;
        }
    }
    for x in &r.edges {
        let e = m.entry((x.row, x.col, x.id)).or_default();
        e.0 += 1;
        e.1 += x.fraction as f64;
    }
    m
}

fn one_record_per_cell(r: &BurnResult) -> bool {
    cells(r).values().all(|&(n, _)| n == 1)
}

fn area(r: &BurnResult) -> f64 {
    cells(r).values().map(|v| v.1).sum()
}

#[test]
fn parts_touching_mid_cell_make_full_cells() {
    // a shared border at x = 5.5 splits column 5 into two halves
    let g = multi(vec![rect(1.0, 1.0, 5.5, 9.0), rect(5.5, 1.0, 9.0, 9.0)]);
    let r = burn(std::slice::from_ref(&g), &g10(), BurnOptions::coverage()).unwrap();
    assert!(one_record_per_cell(&r));
    assert!((area(&r) - 64.0).abs() < 1e-4);
    // column 5 is covered by runs, rows 1..9
    assert!(r.edges.iter().all(|e| e.col != 5));
    // and the runs coalesce across the old border: one run per row
    assert_eq!(r.runs.len(), 8);

    // parity mode keeps the parts' own records, as the C++ core
    let p = burn(std::slice::from_ref(&g), &g10(), BurnOptions::parity(BurnMode::Coverage)).unwrap();
    assert!(!one_record_per_cell(&p));
    assert!((area(&p) - 64.0).abs() < 1e-4);
}

#[test]
fn two_parts_in_one_cell_make_one_edge() {
    let g = multi(vec![rect(4.1, 4.1, 4.4, 4.9), rect(4.6, 4.1, 4.9, 4.9)]);
    let r = burn(&[g], &g10(), BurnOptions::coverage()).unwrap();
    assert_eq!(r.edges.len(), 1);
    assert!(r.runs.is_empty());
    assert!((r.edges[0].fraction as f64 - 0.48).abs() < 1e-6);
}

#[test]
fn overlapping_parts_never_exceed_a_full_cell() {
    // the overlap is not unioned, but no cell holds more than one record
    // or more than a whole cell
    let g = multi(vec![rect(1.0, 1.0, 6.3, 6.3), rect(3.7, 3.7, 8.6, 8.6)]);
    for mode in [BurnMode::Coverage, BurnMode::Approx] {
        let r = burn(std::slice::from_ref(&g), &g10(), BurnOptions::default().with_mode(mode)).unwrap();
        assert!(one_record_per_cell(&r), "{mode:?}");
        assert!(cells(&r).values().all(|&(_, v)| v <= 1.0 + 1e-6), "{mode:?}");
        // the union of the parts is between the larger part and their sum
        let a = area(&r);
        assert!((5.3 * 5.3..=5.3 * 5.3 + 4.9 * 4.9).contains(&a), "{mode:?} {a}");
    }
}

#[test]
fn disjoint_parts_keep_their_area() {
    // the merge only touches cells that several parts cross
    let parts = vec![rect(0.3, 0.3, 3.2, 2.7), rect(3.6, 0.3, 6.1, 4.4), rect(1.1, 5.2, 8.8, 9.6)];
    let want: f64 = 2.9 * 2.4 + 2.5 * 4.1 + 7.7 * 4.4;
    let r = burn(&[multi(parts.clone())], &g10(), BurnOptions::coverage()).unwrap();
    assert!(one_record_per_cell(&r));
    assert!((area(&r) - want).abs() < 1e-4);
    // the same cells, with the same values, as the parts burned one by one
    let singles: Vec<Geometry> = parts.into_iter().map(Geometry::Polygon).collect();
    let s = burn(&singles, &g10(), BurnOptions::coverage()).unwrap();
    let mut by_cell: HashMap<(i32, i32), f64> = HashMap::new();
    for ((row, col, _), (_, v)) in cells(&s) {
        *by_cell.entry((row, col)).or_default() += v;
    }
    let merged: HashMap<(i32, i32), f64> =
        cells(&r).into_iter().map(|((row, col, _), (_, v))| ((row, col), v)).collect();
    assert_eq!(merged.len(), by_cell.len());
    for (k, v) in by_cell {
        assert!((merged[&k] - v).abs() < 1e-6, "{k:?}");
    }
}

#[test]
fn line_parts_sum_per_cell() {
    let part = ring(&[(0.5, 2.5), (9.5, 2.5)]);
    let g = Geometry::MultiLineString(vec![part.clone(), part]);
    let r = burn(std::slice::from_ref(&g), &g10(), BurnOptions::default()).unwrap();
    assert_eq!(r.lines.len(), 10);
    let total: f64 = r.lines.iter().map(|l| l.length as f64).sum();
    assert!((total - 18.0).abs() < 1e-5);
    let p = burn(std::slice::from_ref(&g), &g10(), BurnOptions::parity(BurnMode::Coverage)).unwrap();
    assert_eq!(p.lines.len(), 20);
}

#[test]
fn short_rings_and_lines_are_noted() {
    let good = Geometry::Polygon(rect(1.0, 1.0, 3.0, 3.0));
    let short_hole = Geometry::Polygon(Polygon::new(vec![
        rect(1.0, 1.0, 9.0, 9.0).rings[0].clone(),
        ring(&[(4.0, 4.0), (5.0, 5.0)]),
    ]));
    let short_part = multi(vec![Polygon::new(vec![ring(&[(1.0, 1.0), (2.0, 2.0)])]), rect(5.0, 5.0, 7.0, 7.0)]);
    let short_line = Geometry::MultiLineString(vec![ring(&[(1.0, 1.0)]), ring(&[(0.5, 0.5), (9.5, 0.5)])]);
    let r = burn(&[good.clone(), short_hole, short_part, short_line, good], &g10(), BurnOptions::default()).unwrap();
    let notes: Vec<(i32, &str)> = r.notes.iter().map(|n| (n.geom_index, n.message.as_str())).collect();
    assert_eq!(notes.len(), 3, "{notes:?}");
    assert_eq!(notes[0].0, 1);
    assert!(notes[0].1.contains("hole 1 of polygon part 0"), "{}", notes[0].1);
    assert_eq!(notes[1].0, 2);
    assert!(notes[1].1.contains("polygon part 0: exterior ring"), "{}", notes[1].1);
    assert_eq!(notes[2].0, 3);
    assert!(notes[2].1.contains("line part 0"), "{}", notes[2].1);
    // the rest of each geometry is burned: the exterior without its hole,
    // the other polygon part, the other line part, and both good squares
    let by_id = |id: i32| cells(&r).iter().filter(|(k, _)| k.2 == id).map(|(_, v)| v.1).sum::<f64>();
    assert!((by_id(0) - 4.0).abs() < 1e-5 && (by_id(4) - 4.0).abs() < 1e-5);
    assert!((by_id(1) - 64.0).abs() < 1e-5);
    assert!((by_id(2) - 4.0).abs() < 1e-5);
    assert_eq!(r.lines.iter().filter(|l| l.id == 3).count(), 10);
}

#[test]
fn point_empty_is_empty() {
    // WKB POINT EMPTY: little-endian point, both coordinates NaN
    let mut empty = vec![1u8, 1, 0, 0, 0];
    empty.extend_from_slice(&f64::NAN.to_le_bytes());
    empty.extend_from_slice(&f64::NAN.to_le_bytes());
    let r = burn_wkb([empty.as_slice()], &g10(), BurnOptions::default()).unwrap();
    assert!(r.is_empty() && r.notes.is_empty(), "{:?}", r.notes);
    let nan_x = Geometry::Point(Coord::new(f64::NAN, 1.0));
    let r = burn(&[nan_x], &g10(), BurnOptions::default()).unwrap();
    assert_eq!(r.notes.len(), 1, "a point with one NaN coordinate is not empty");
}

#[test]
fn id_offset_concatenates() {
    let geoms: Vec<Geometry> = (0..6)
        .map(|i| Geometry::Polygon(rect(0.3 + i as f64, 0.3 + i as f64, 3.6 + i as f64, 2.2 + i as f64)))
        .collect();
    let whole = burn(&geoms, &g10(), BurnOptions::default()).unwrap();
    let a = burn(&geoms[..4], &g10(), BurnOptions::default()).unwrap();
    let b = burn(&geoms[4..], &g10(), BurnOptions::default().with_id_offset(4)).unwrap();
    let mut joined = a;
    joined.runs.extend(b.runs);
    joined.edges.extend(b.edges);
    assert_eq!(joined, whole);

    // notes carry the offset id too
    let short = Geometry::Polygon(Polygon::new(vec![ring(&[(1.0, 1.0), (2.0, 2.0)])]));
    let r = burn(&[short], &g10(), BurnOptions::default().with_id_offset(7)).unwrap();
    assert_eq!(r.notes[0].geom_index, 7);
}

#[test]
fn id_offset_out_of_range() {
    let g = Geometry::Polygon(rect(1.0, 1.0, 3.0, 3.0));
    let neg = burn(std::slice::from_ref(&g), &g10(), BurnOptions::default().with_id_offset(-1));
    assert!(matches!(neg, Err(BurnError::InvalidOptions(_))));
    let two = [g.clone(), g];
    assert!(burn(&two[..1], &g10(), BurnOptions::default().with_id_offset(i32::MAX)).is_ok());
    let over = burn(&two, &g10(), BurnOptions::default().with_id_offset(i32::MAX));
    assert!(matches!(over, Err(BurnError::InvalidOptions(_))));
}
