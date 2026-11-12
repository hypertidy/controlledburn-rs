//! One record per cell and id, for a whole result.

use crate::coalesce::coalesce_runs;
use crate::output::{BurnResult, GridRun};

/// The sweep's tolerance on a cell's coverage, as in `burn`.
const TOL: f32 = 1e-6;

/// Make `result` hold one record per `(id, row, col)`, by the rule
/// [`burn`](crate::burn) uses for the parts of a multi-part geometry
/// ([`BurnOptions::merge_parts`](crate::BurnOptions)).
///
/// - Edges in the same cell with the same id become one edge, their
///   fractions summed (`f32`, in table order), or a one-cell run once
///   the sum reaches `1 - 1e-6`.
/// - Runs are sorted and coalesced per id and row
///   ([`coalesce_runs`]); overlapping runs become one.
/// - An edge in a cell a run of the same id covers is dropped.
/// - Lines in the same cell with the same id become one line, their
///   lengths summed in table order.
///
/// Points and notes are left as they are. Overlaps are not unioned: a
/// cell's sum is capped at a full cell, as in `burn`.
///
/// This is for records that come together after the burn, such as a
/// front end that burns in an unwrapped index space and folds columns
/// or rows back (a periodic or tripolar grid), where two pieces of one
/// geometry can land in the same cell. For a single burn with
/// `merge_parts` on (the default) it changes nothing but the order.
///
/// The tables come out sorted by `(id, row, col)` (runs by `col_start`).
///
/// ```
/// use controlledburn::{merge_cells, BurnResult, GridEdge};
///
/// // one geometry, two pieces folded onto the same cell
/// let mut r = BurnResult::default();
/// r.edges.push(GridEdge { row: 2, col: 0, fraction: 0.25, id: 7 });
/// r.edges.push(GridEdge { row: 2, col: 0, fraction: 0.5, id: 7 });
/// merge_cells(&mut r);
/// assert_eq!(r.edges, vec![GridEdge { row: 2, col: 0, fraction: 0.75, id: 7 }]);
/// ```
pub fn merge_cells(result: &mut BurnResult) {
    // stable: a cell's fractions are summed in table order
    let edges = &mut result.edges;
    edges.sort_by_key(|e| (e.id, e.row, e.col));
    let mut w = 0;
    for i in 0..edges.len() {
        let cur = edges[i];
        if w > 0 {
            let last = &mut edges[w - 1];
            if (last.id, last.row, last.col) == (cur.id, cur.row, cur.col) {
                last.fraction += cur.fraction;
                continue;
            }
        }
        edges[w] = cur;
        w += 1;
    }
    edges.truncate(w);

    let mut w = 0;
    for i in 0..edges.len() {
        let e = edges[i];
        if e.fraction >= 1.0 - TOL {
            result.runs.push(GridRun { row: e.row, col_start: e.col, col_end: e.col + 1, id: e.id });
        } else {
            edges[w] = e;
            w += 1;
        }
    }
    edges.truncate(w);
    coalesce_runs(&mut result.runs);

    // both sorted by (id, row, col): drop the edges a run covers
    let runs = &result.runs;
    let mut j = 0;
    let mut w = 0;
    for i in 0..edges.len() {
        let e = edges[i];
        while j < runs.len() && (runs[j].id, runs[j].row, runs[j].col_end) <= (e.id, e.row, e.col) {
            j += 1;
        }
        let covered = j < runs.len() && (runs[j].id, runs[j].row) == (e.id, e.row) && runs[j].col_start <= e.col;
        if !covered {
            edges[w] = e;
            w += 1;
        }
    }
    edges.truncate(w);

    let lines = &mut result.lines;
    lines.sort_by_key(|l| (l.id, l.row, l.col));
    let mut w = 0;
    for i in 0..lines.len() {
        let cur = lines[i];
        if w > 0 {
            let last = &mut lines[w - 1];
            if (last.id, last.row, last.col) == (cur.id, cur.row, cur.col) {
                last.length += cur.length;
                continue;
            }
        }
        lines[w] = cur;
        w += 1;
    }
    lines.truncate(w);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::output::{GridEdge, GridLine};

    fn e(row: i32, col: i32, fraction: f32, id: i32) -> GridEdge {
        GridEdge { row, col, fraction, id }
    }

    #[test]
    fn a_sum_within_tolerance_of_one_is_a_full_cell() {
        let mut r = BurnResult {
            edges: vec![e(0, 3, 0.5, 1), e(0, 3, 0.4999995, 1), e(0, 4, 0.5, 1), e(0, 4, 0.499, 1)],
            ..Default::default()
        };
        merge_cells(&mut r);
        assert_eq!(r.runs, vec![GridRun { row: 0, col_start: 3, col_end: 4, id: 1 }]);
        assert_eq!(r.edges.len(), 1);
        assert!((r.edges[0].fraction - 0.999).abs() < 1e-6);
    }

    #[test]
    fn folded_pieces_make_one_record_per_cell_and_id() {
        // a geometry's run and edges folded onto cells it already has,
        // beside another id's records in the same cells
        let mut r = BurnResult {
            runs: vec![
                GridRun { row: 1, col_start: 0, col_end: 3, id: 0 },
                GridRun { row: 1, col_start: 2, col_end: 5, id: 0 },
            ],
            edges: vec![e(1, 4, 0.3, 0), e(1, 5, 0.3, 0), e(1, 5, 0.2, 0), e(1, 5, 0.6, 1)],
            lines: vec![
                GridLine { row: 2, col: 2, length: 0.5, id: 3 },
                GridLine { row: 2, col: 2, length: 0.25, id: 3 },
                GridLine { row: 2, col: 2, length: 1.0, id: 4 },
            ],
            ..Default::default()
        };
        merge_cells(&mut r);
        assert_eq!(r.runs, vec![GridRun { row: 1, col_start: 0, col_end: 5, id: 0 }]);
        assert_eq!(r.edges, vec![e(1, 5, 0.5, 0), e(1, 5, 0.6, 1)]);
        let lines: Vec<(i32, f32)> = r.lines.iter().map(|l| (l.id, l.length)).collect();
        assert_eq!(lines, vec![(3, 0.75), (4, 1.0)]);
    }
}
