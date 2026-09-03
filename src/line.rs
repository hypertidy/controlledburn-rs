//! Line rasterization: length of the line within each cell it touches.

use crate::geometry::Coord;
use crate::grid::Grid;
use crate::output::GridLine;
use crate::walker::walk_polyline;

/// Walk the line on the full padded grid (no sub-grid: axis-parallel
/// lines have a degenerate envelope, and the walker's cost is O(cells
/// touched) regardless of grid size) and emit one record per in-grid
/// cell with the summed segment length of every traversal through it.
pub(crate) fn process_line(
    line: &[Coord],
    full_grid: &Grid<0>,
    line_id: i32,
    out: &mut Vec<GridLine>,
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

    let cells = walk_polyline(line.to_vec(), &grid, false)?;

    let first = out.len();
    let arena = cells.coords.as_slice();
    for cr in &cells.records {
        let (r, c) = (cr.row, cr.col);
        // Skip the padding ring.
        if r < 1 || r >= n_rows - 1 || c < 1 || c >= n_cols - 1 {
            continue;
        }
        let mut total_length = 0.0f64;
        for t in &cr.traversals {
            let tc = t.coords(arena);
            if tc.len() < 2 {
                continue;
            }
            for w in tc.windows(2) {
                let sx = w[1].x - w[0].x;
                let sy = w[1].y - w[0].y;
                total_length += (sx * sx + sy * sy).sqrt();
            }
        }
        if total_length > 0.0 {
            // Padded (r, c) map to 1-based bounded indices directly:
            // padding = 1 absorbs the +1.
            out.push(GridLine { row: r as i32, col: c as i32, length: total_length as f32, id: line_id });
        }
    }
    // Row-major record order, as the C++ (which iterates an ordered map).
    out[first..].sort_by_key(|l| (l.row, l.col));
    Ok(())
}
