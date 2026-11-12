//! Line rasterization: length of the line within each cell it touches.

use crate::ee::bbox::BBox;
use crate::geometry::Coord;
use crate::grid::Grid;
use crate::output::GridLine;
use crate::walker::walk_polyline;

/// Walk the line on the full padded grid (no sub-grid: axis-parallel
/// lines have a degenerate envelope, and the walker's cost is O(cells
/// touched) regardless of grid size) and emit one record per in-grid
/// cell with the summed segment length of every traversal through it.
///
/// `near` (a window widened by two cells): walk only the parts of the
/// line near it ([`near_runs`]); a cell's lengths are summed over the
/// parts in line order, so the cells near the window get the whole walk's
/// lengths, bit for bit.
pub(crate) fn process_line(
    line: &[Coord],
    full_grid: &Grid<0>,
    line_id: i32,
    out: &mut Vec<GridLine>,
    near: Option<&BBox>,
) -> Result<(), String> {
    if line.len() < 2 {
        return Ok(());
    }
    let grid = full_grid.make_infinite();
    if grid.is_empty() {
        return Ok(());
    }
    let n_rows = grid.rows();
    let n_cols = grid.cols();

    let runs: Vec<&[Coord]> = match near {
        None => vec![line],
        Some(b) => near_runs(line, full_grid, b),
    };

    // (padded row, col) -> summed length, in walk order. One walk visits
    // a cell in one record; several runs can each visit it.
    let mut lengths: Vec<((usize, usize), f64)> = Vec::new();
    let mut index = std::collections::HashMap::new();
    let several = runs.len() > 1;
    for run in runs {
        let cells = walk_polyline(run.to_vec(), &grid, false)?;
        let arena = cells.coords.as_slice();
        for cr in &cells.records {
            let (r, c) = (cr.row, cr.col);
            // Skip the padding ring.
            if r < 1 || r >= n_rows - 1 || c < 1 || c >= n_cols - 1 {
                continue;
            }
            let k = if several {
                *index.entry((r, c)).or_insert_with(|| {
                    lengths.push(((r, c), 0.0));
                    lengths.len() - 1
                })
            } else {
                lengths.push(((r, c), 0.0));
                lengths.len() - 1
            };
            let total_length = &mut lengths[k].1;
            for t in cells.traversals(cr) {
                let tc = t.coords(arena);
                if tc.len() < 2 {
                    continue;
                }
                for w in tc.windows(2) {
                    let sx = w[1].x - w[0].x;
                    let sy = w[1].y - w[0].y;
                    *total_length += (sx * sx + sy * sy).sqrt();
                }
            }
        }
    }
    let first = out.len();
    for ((r, c), total_length) in lengths {
        if total_length > 0.0 {
            // Strip the padding ring to get 0-based grid indices.
            out.push(GridLine { row: (r - 1) as i32, col: (c - 1) as i32, length: total_length as f32, id: line_id });
        }
    }
    // Row-major record order, as the C++ (which iterates an ordered map).
    out[first..].sort_by_key(|l| (l.row, l.col));
    Ok(())
}

/// The parts of `line` to walk for the cells in `near`: each run of
/// segments whose box meets `near`, started at the last vertex before it
/// that lies inside one cell, clear of every wall (or at the line's
/// start), and runs that then overlap joined.
///
/// From such a vertex on, a walk started there and the whole walk are the
/// same: they are in the same cell, and the walker computes each crossing
/// from the segment's own vertices, not from where it entered the cell.
/// From a vertex on a wall they need not be: a segment along a wall stays
/// on the side the walk arrived from, so a walk started at a corner can
/// put a wall segment in the cell across it. The vertex the walk starts
/// at, and the one it ends at, are outside `near` (they belong to a
/// segment that misses it), so the cells where the part's walk differs
/// from the whole walk are not the window's.
fn near_runs<'a>(line: &'a [Coord], grid: &Grid<0>, b: &BBox) -> Vec<&'a [Coord]> {
    let meets = |p: &Coord, q: &Coord| {
        !(p.x.max(q.x) < b.xmin || p.x.min(q.x) > b.xmax || p.y.max(q.y) < b.ymin || p.y.min(q.y) > b.ymax)
    };
    let (dx, dy) = (grid.dx(), grid.dy());
    // well clear of the walls: of rounding in `xmin + k * dx` and of the
    // walker's snap (1e-9 of a cell)
    let clear = |v: f64, d: f64| {
        let f = v / d;
        (f - f.round()).abs() * d > 1e-6 * d + 1e-13 * v.abs()
    };
    let inside_one_cell = |p: &Coord| clear(p.x - grid.xmin(), dx) && clear(grid.ymax() - p.y, dy);
    let mut runs: Vec<(usize, usize)> = Vec::new();
    let mut i = 0;
    let n = line.len() - 1;
    while i < n {
        if !meets(&line[i], &line[i + 1]) {
            i += 1;
            continue;
        }
        let mut start = i;
        while start > 0 && !inside_one_cell(&line[start]) {
            start -= 1;
        }
        let mut end = i + 1;
        while end < n && meets(&line[end], &line[end + 1]) {
            end += 1;
        }
        match runs.last_mut() {
            Some(last) if start <= last.1 => last.1 = end,
            _ => runs.push((start, end)),
        }
        i = end;
    }
    runs.into_iter().map(|(s, e)| &line[s..=e]).collect()
}
