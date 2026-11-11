//! An external oracle: burnscars' 50 fixture cases, with coverage from
//! exactextract 0.3.0 and the cell-centre rule from GDAL's rasterize
//! (`fixtures/oracle/`, see its README).
//!
//! Per case, the cells where cb and the oracle disagree (a cell on one
//! side only, or a fraction more than 1e-6 apart) must be exactly the
//! expected set: the cells burnscars lists in `known_diffs.toml`, less the
//! approx cells where cb agrees with GDAL and burnscars does not
//! (`CB_AGREES_WITH_GDAL`). A new difference fails, and so does a listed
//! one that goes away, so the list stays true.

use std::collections::{BTreeMap, BTreeSet};

use controlledburn::{burn_wkb, BurnMode, BurnOptions, BurnResult, GridSpec};

const CASES: &str = include_str!("../fixtures/oracle/cases.csv");
const COVERAGE: &str = include_str!("../fixtures/oracle/coverage.csv");
const APPROX: &str = include_str!("../fixtures/oracle/approx.csv");
const KNOWN: &str = include_str!("../fixtures/oracle/known_diffs.toml");
const TOL: f64 = 1e-6;

type Cell = (u32, u32);

/// Approx cells burnscars lists as differing from GDAL where cb gives
/// GDAL's answer. All are cell centres within about 4e-14 of an edge:
/// cb's exact side test (`orient2d`) decides them as GDAL does, and
/// burnscars' floating-point intercept does not.
const CB_AGREES_WITH_GDAL: &[(&str, &[Cell])] = &[
    ("random_9", &[(7, 17)]),
    ("random_24", &[(12, 19), (13, 20), (20, 34)]),
    ("random_27", &[(5, 5)]),
    ("random_30", &[(22, 5)]),
    ("random_39", &[(37, 15), (46, 22)]),
];

fn rows(csv: &str) -> Vec<Vec<&str>> {
    csv.lines().skip(1).filter(|l| !l.is_empty()).map(|l| l.split(',').collect()).collect()
}

fn hex(s: &str) -> Vec<u8> {
    (0..s.len()).step_by(2).map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap()).collect()
}

/// `known_diffs.toml`: (table, case) -> cells. Only the shape the file
/// uses: `[[diff]]` blocks of `key = value` lines, `cells` on one line.
fn known_diffs() -> BTreeMap<(String, String), BTreeSet<Cell>> {
    let mut out = BTreeMap::new();
    let (mut table, mut case) = (String::new(), String::new());
    let unquote = |v: &str| v.trim().trim_matches('"').to_string();
    for line in KNOWN.lines() {
        let Some((k, v)) = line.split_once(" = ") else { continue };
        match k.trim() {
            "table" => table = unquote(v),
            "case" => case = unquote(v),
            "cells" => {
                let nums: Vec<u32> = v
                    .split(|c: char| !c.is_ascii_digit())
                    .filter(|s| !s.is_empty())
                    .map(|s| s.parse().unwrap())
                    .collect();
                let cells = nums.chunks(2).map(|p| (p[0], p[1])).collect();
                out.insert((table.clone(), case.clone()), cells);
            }
            _ => {}
        }
    }
    out
}

fn burn_case(case: &[&str], mode: BurnMode) -> BurnResult {
    let f = |i: usize| case[i].parse::<f64>().unwrap();
    let grid = GridSpec::new(f(2), f(3), f(4), f(5), case[6].parse().unwrap(), case[7].parse().unwrap());
    let r = burn_wkb([hex(case[1]).as_slice()], &grid, BurnOptions::default().with_mode(mode)).unwrap();
    assert!(r.notes.is_empty(), "{}: {:?}", case[0], r.notes);
    r
}

/// Per-cell value (runs 1, edges their fraction, summed), cells below
/// TOL dropped.
fn cells(r: &BurnResult) -> BTreeMap<Cell, f64> {
    let mut m: BTreeMap<Cell, f64> = BTreeMap::new();
    for x in &r.runs {
        for c in x.col_start..x.col_end {
            *m.entry((x.row as u32, c as u32)).or_default() += 1.0;
        }
    }
    for e in &r.edges {
        *m.entry((e.row as u32, e.col as u32)).or_default() += e.fraction as f64;
    }
    m.retain(|_, f| *f >= TOL);
    m
}

fn oracle_cells(rows: &[Vec<&str>], case: &str) -> BTreeMap<Cell, f64> {
    rows.iter()
        .filter(|r| r[0] == case)
        .map(|r| ((r[1].parse().unwrap(), r[2].parse().unwrap()), r.get(3).map_or(1.0, |f| f.parse().unwrap())))
        .filter(|(_, f)| *f >= TOL)
        .collect()
}

fn diff_cells(got: &BTreeMap<Cell, f64>, want: &BTreeMap<Cell, f64>) -> BTreeSet<Cell> {
    got.keys()
        .chain(want.keys())
        .filter(|k| match (got.get(k), want.get(k)) {
            (Some(a), Some(b)) => (a - b).abs() > TOL,
            _ => true,
        })
        .copied()
        .collect()
}

fn check(table: &str, mode: BurnMode, oracle: &str) {
    let known = known_diffs();
    let oracle = rows(oracle);
    let cases = rows(CASES);
    assert_eq!(cases.len(), 50);
    let mut failures = Vec::new();
    for case in &cases {
        let name = case[0];
        let got = diff_cells(&cells(&burn_case(case, mode)), &oracle_cells(&oracle, name));
        let mut want = known.get(&(table.to_string(), name.to_string())).cloned().unwrap_or_default();
        if mode == BurnMode::Approx {
            for (_, agree) in CB_AGREES_WITH_GDAL.iter().filter(|(c, _)| *c == name) {
                for cell in *agree {
                    // listed by burnscars, so it must be there to remove
                    assert!(want.remove(cell), "{name} {cell:?} is not in known_diffs.toml");
                }
            }
        }
        if got != want {
            failures.push(format!("{table} {name}: differs at {got:?}, expected {want:?}"));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn coverage_matches_exactextract() {
    check("coverage", BurnMode::Coverage, COVERAGE);
}

#[test]
fn approx_matches_gdal_centre_rule() {
    check("approx", BurnMode::Approx, APPROX);
}
