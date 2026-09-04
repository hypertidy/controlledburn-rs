//! Run coalescing: merge adjacent runs on the same row with the same id.
//!
//! The C++ sweep emits one `GridRun` per fully covered boundary cell and
//! one per interior span between boundary cells, without merging
//! neighbours. This pass merges them. It is applied to a polygon's runs
//! after its sweep (runs from one polygon are emitted row by row, left to
//! right) so a linear scan suffices; the general entry point sorts first.

use crate::output::GridRun;

/// Merge adjacent runs in place, assuming `runs` is already ordered by
/// (id, row, col_start) with non-overlapping ranges, as the sweep emits
/// them for a single polygon.
pub(crate) fn coalesce_sorted(runs: &mut Vec<GridRun>) {
    let mut w = 0usize;
    for i in 0..runs.len() {
        let cur = runs[i];
        if w > 0 {
            let last = &mut runs[w - 1];
            if last.id == cur.id && last.row == cur.row && cur.col_start <= last.col_end {
                if cur.col_end > last.col_end {
                    last.col_end = cur.col_end;
                }
                continue;
            }
        }
        runs[w] = cur;
        w += 1;
    }
    runs.truncate(w);
}

/// Merge adjacent or overlapping runs with the same id on the same row.
/// The output is sorted by (id, row, col_start).
pub fn coalesce_runs(runs: &mut Vec<GridRun>) {
    runs.sort_by_key(|r| (r.id, r.row, r.col_start));
    coalesce_sorted(runs);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(row: i32, a: i32, b: i32, id: i32) -> GridRun {
        GridRun { row, col_start: a, col_end: b, id }
    }

    #[test]
    fn merges_neighbours_only() {
        let mut v = vec![r(1, 3, 4, 1), r(1, 4, 5, 1), r(1, 5, 6, 1), r(1, 7, 10, 1), r(2, 1, 2, 1), r(1, 6, 7, 2)];
        coalesce_runs(&mut v);
        assert_eq!(v, vec![r(1, 3, 6, 1), r(1, 7, 10, 1), r(2, 1, 2, 1), r(1, 6, 7, 2)]);
    }
}
