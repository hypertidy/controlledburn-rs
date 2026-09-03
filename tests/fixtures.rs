//! The cross-language aggregate parity contract from
//! fixtures/{geometries,expected}.csv (shared with the C++, R and Python
//! surfaces of the original project).

mod common;

use common::*;
use controlledburn::{burn_wkb, BurnOptions};
use std::collections::HashMap;

fn split_csv(line: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut q = false;
    for ch in line.chars() {
        match ch {
            '"' => q = !q,
            ',' if !q => {
                out.push(std::mem::take(&mut cur));
            }
            _ => cur.push(ch),
        }
    }
    out.push(cur);
    out
}

fn read_csv(name: &str) -> Vec<HashMap<String, String>> {
    let text = std::fs::read_to_string(fixtures_dir().join(name)).unwrap();
    let mut lines = text.lines();
    let header = split_csv(lines.next().unwrap());
    lines.filter(|l| !l.trim().is_empty()).map(|l| header.iter().cloned().zip(split_csv(l)).collect()).collect()
}

#[test]
fn parity_fixtures() {
    let geoms = read_csv("geometries.csv");
    let expected: HashMap<String, HashMap<String, String>> =
        read_csv("expected.csv").into_iter().map(|r| (r["case"].clone(), r)).collect();
    assert_eq!(geoms.len(), 10);

    for g in &geoms {
        let case = &g["case"];
        let e = &expected[case];
        let grid = controlledburn::GridSpec::new(
            g["xmin"].parse().unwrap(),
            g["ymin"].parse().unwrap(),
            g["xmax"].parse().unwrap(),
            g["ymax"].parse().unwrap(),
            g["ncol"].parse().unwrap(),
            g["nrow"].parse().unwrap(),
        );
        let wkb = hex(&g["wkb_hex"]);
        let r = burn_wkb([wkb.as_slice()], &grid, BurnOptions::default()).unwrap();
        assert!(r.notes.is_empty(), "{case}: {:?}", r.notes);

        let near = |got: f64, want: f64, tol: f64| {
            if want != 0.0 {
                ((got - want) / want).abs() <= tol
            } else {
                (got - want).abs() <= tol
            }
        };

        if let Ok(want) = e["covered_area"].parse::<f64>() {
            let tol: f64 = e["tol_rel"].parse().unwrap();
            let got = area(&r, &grid);
            assert!(near(got, want, tol), "{case}: covered area {got} vs {want}");
            // 1 = edges must be empty; 0 = not checked (same as the C++ test).
            if e["edges_empty"] == "1" {
                assert!(r.edges.is_empty(), "{case}: edges must be empty");
            }
        }
        if let Ok(want) = e["line_length"].parse::<f64>() {
            let tol: f64 = e["tol_rel"].parse().unwrap();
            let got = r.line_length();
            assert!(near(got, want, tol), "{case}: line length {got} vs {want}");
        }
        if let Ok(want) = e["n_points"].parse::<usize>() {
            assert_eq!(r.points.len(), want, "{case}: n_points");
        }
    }
}
