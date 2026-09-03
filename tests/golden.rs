//! Exact-table parity against the C++ core's output (both modes) for the
//! parity fixtures and the edge-case probes, in parity mode (coalescing
//! off). Then: applying the coalesce pass to the C++ tables must equal a
//! default-mode run.

mod common;

use common::*;
use controlledburn::{burn, coalesce_runs, BurnMode, BurnOptions, BurnResult};

fn assert_tables_equal(case: &str, mode: &str, got: &BurnResult, want: &BurnResult) {
    assert_eq!(got.runs, want.runs, "{case} [{mode}] runs differ");
    assert_eq!(got.edges.len(), want.edges.len(), "{case} [{mode}] edge count differs");
    for (i, (g, w)) in got.edges.iter().zip(&want.edges).enumerate() {
        assert_eq!((g.row, g.col, g.id), (w.row, w.col, w.id), "{case} [{mode}] edge {i} position differs");
        assert_eq!(g.fraction.to_bits(), w.fraction.to_bits(), "{case} [{mode}] edge {i} fraction differs: {} vs {}", g.fraction, w.fraction);
    }
    assert_eq!(got.lines.len(), want.lines.len(), "{case} [{mode}] line count differs");
    for (i, (g, w)) in got.lines.iter().zip(&want.lines).enumerate() {
        assert_eq!((g.row, g.col, g.id), (w.row, w.col, w.id), "{case} [{mode}] line {i} position differs");
        assert_eq!(g.length.to_bits(), w.length.to_bits(), "{case} [{mode}] line {i} length differs: {} vs {}", g.length, w.length);
    }
    assert_eq!(got.points, want.points, "{case} [{mode}] points differ");
    assert!(got.notes.is_empty(), "{case} [{mode}] unexpected notes: {:?}", got.notes);
}

#[test]
fn golden_parity_both_modes() {
    let cases = load_golden();
    assert_eq!(cases.len(), 25);
    for c in &cases {
        let name = c["case"].as_str().unwrap();
        let g = golden_geometry(&c["geometry"]);
        let grid = golden_grid(c);
        for (mode_name, mode) in [("coverage", BurnMode::Coverage), ("approx", BurnMode::Approx)] {
            let got = burn(std::slice::from_ref(&g), &grid, BurnOptions::parity(mode)).unwrap();
            let want = golden_result(&c[mode_name]);
            assert_tables_equal(name, mode_name, &got, &want);
        }
    }
}

#[test]
fn golden_parity_via_wkb() {
    for c in &load_golden() {
        let Some(h) = c["wkb_hex"].as_str().filter(|s| !s.is_empty()) else { continue };
        let bytes = hex(h);
        let grid = golden_grid(c);
        let got = controlledburn::burn_wkb([bytes.as_slice()], &grid, BurnOptions::parity(BurnMode::Coverage)).unwrap();
        assert_tables_equal(c["case"].as_str().unwrap(), "coverage/wkb", &got, &golden_result(&c["coverage"]));
    }
}

#[test]
fn coalescing_matches_post_pass_on_golden() {
    for c in &load_golden() {
        let name = c["case"].as_str().unwrap();
        let g = golden_geometry(&c["geometry"]);
        let grid = golden_grid(c);
        for (mode_name, mode) in [("coverage", BurnMode::Coverage), ("approx", BurnMode::Approx)] {
            let mut want = golden_result(&c[mode_name]);
            coalesce_runs(&mut want.runs);
            let got = burn(std::slice::from_ref(&g), &grid, BurnOptions::default().with_mode(mode)).unwrap();
            let mut got_runs = got.runs.clone();
            got_runs.sort_by_key(|r| (r.id, r.row, r.col_start));
            assert_eq!(got_runs, want.runs, "{name} [{mode_name}] coalesced runs differ");
            assert_eq!(got.edges.len(), want.edges.len());
            assert!(got.runs.len() <= want.runs.len() || got.runs == want.runs);
        }
    }
}
