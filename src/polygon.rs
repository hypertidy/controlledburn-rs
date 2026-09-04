//! Coverage-mode polygon rasterization: walk each ring, accumulate
//! per-boundary-cell coverage fractions and winding deltas, then sweep
//! each row to emit interior runs and boundary edges.
//!
//! Each polygon (each MultiPolygon component) is processed independently
//! with its own sub-grid and winding sweep, so winding from one disjoint
//! component never bleeds into another's boundary cells.

use crate::coverage::{analytical_covered_fraction, closed_ring_covered_fraction};
use crate::ee::bbox::BBox;
use crate::ee::traversal_areas::{left_hand_area_with, Scratch};
use crate::geometry::{is_ccw, Coord, Polygon};
use crate::grid::Grid;
use crate::output::{GridEdge, GridRun};
use crate::walker::{walk_polyline, Traversal};

/// Per-cell boundary data for the winding sweep. `col` is the 0-based
/// column in the full grid; padding columns are -1 and ncol.
#[derive(Clone, Copy, Debug)]
struct BoundaryCellRecord {
    col: i64,
    /// Accumulated (signed) coverage fraction. Kept in f32 to match the
    /// C++ core bit for bit.
    coverage: f32,
    winding_delta: i32,
}

pub(crate) fn ring_envelope(ring: &[Coord]) -> BBox {
    if ring.is_empty() {
        return BBox::make_empty();
    }
    let (mut xmin, mut xmax, mut ymin, mut ymax) = (ring[0].x, ring[0].x, ring[0].y, ring[0].y);
    for c in ring {
        if c.x < xmin {
            xmin = c.x;
        }
        if c.x > xmax {
            xmax = c.x;
        }
        if c.y < ymin {
            ymin = c.y;
        }
        if c.y > ymax {
            ymax = c.y;
        }
    }
    BBox::new(xmin, ymin, xmax, ymax)
}

/// Walk one ring and add its boundary-cell coverage and winding deltas
/// to `row_data` (indexed by sub-grid row).
#[allow(clippy::too_many_arguments)]
fn walk_ring(
    mut coords: Vec<Coord>,
    ccw: bool,
    is_exterior: bool,
    grid: &Grid<1>,
    row_data: &mut [Vec<BoundaryCellRecord>],
    sub_rows: usize,
    sub_cols: usize,
    col_off: usize,
    scratch: &mut Scratch,
) -> Result<(), String> {
    if coords.len() < 4 {
        return Ok(());
    }
    // Normalise to CCW so "left of travel" is the interior.
    if !ccw {
        coords.reverse();
    }
    let coverage_factor: f32 = if is_exterior { 1.0 } else { -1.0 };
    let winding_factor: i32 = if is_exterior { 1 } else { -1 };

    let cells = walk_polyline(coords, grid, true)?;
    let arena = cells.coords.as_slice();
    let mut valid: Vec<&Traversal> = Vec::with_capacity(8);
    let mut lists: Vec<&[Coord]> = Vec::with_capacity(8);

    for cr in &cells.records {
        let (r, c) = (cr.row, cr.col);
        // Skip padding ROWS: they affect no grid row's winding.
        if r < 1 {
            continue;
        }
        let sub_r = r - 1;
        if sub_r >= sub_rows {
            continue;
        }
        // Padding COLUMNS still carry winding deltas for their row (an
        // edge entirely outside the grid but crossing rows).
        let (full_col, in_grid_cols): (i64, bool) = if c < 1 {
            (col_off as i64 - 1, false)
        } else {
            let sub_c = c - 1;
            if sub_c >= sub_cols {
                ((col_off + sub_cols) as i64, false)
            } else {
                ((col_off + sub_c) as i64, true)
            }
        };

        // Valid traversals: proper enter + exit with movement, or a closed
        // ring entirely within this cell.
        valid.clear();
        valid.extend(cr.traversals.iter().filter(|t| {
            (t.traversed() && t.multiple_unique_coordinates(arena)) || (t.entry.is_none() && t.is_closed_ring(arena))
        }));
        if valid.is_empty() {
            continue;
        }

        // Coverage fraction, in-grid cells only.
        let mut frac: f32 = 0.0;
        if in_grid_cols {
            if valid.len() == 1 && valid[0].entry.is_none() && valid[0].is_closed_ring(arena) {
                frac = closed_ring_covered_fraction(&cr.bbox, valid[0].coords(arena)) as f32;
            } else if valid.len() == 1 {
                frac = analytical_covered_fraction(&cr.bbox, valid[0].coords(arena)) as f32;
            } else {
                lists.clear();
                lists.extend(valid.iter().map(|t| t.coords(arena)));
                let cell_area = cr.bbox.area();
                if cell_area > 0.0 {
                    frac = (left_hand_area_with(&cr.bbox, &lists, scratch) / cell_area) as f32;
                }
            }
        }

        // Records are pushed unconditionally and merged in the sweep after a
        // stable sort by column. Within one ring a cell is visited once, so
        // the merge adds the same contributions in the same order as the
        // C++ find_or_create accumulation (which was a linear scan per
        // row, quadratic in boundary cells).
        let row_vec = &mut row_data[sub_r];

        if frac != 0.0 {
            row_vec.push(BoundaryCellRecord { col: full_col, coverage: coverage_factor * frac, winding_delta: 0 });
        }

        // Winding deltas must be stored even when coverage is zero: a
        // traversal along a cell wall has zero area but still crosses the
        // row centre.
        let y_mid = (cr.bbox.ymin + cr.bbox.ymax) / 2.0;
        let mut winding_delta = 0;
        for t in &valid {
            let tc = t.coords(arena);
            if !t.traversed() || tc.len() < 2 {
                continue;
            }
            let entry_y = tc[0].y;
            let exit_y = tc[tc.len() - 1].y;
            let crosses = (entry_y > y_mid && exit_y < y_mid) || (entry_y < y_mid && exit_y > y_mid);
            if !crosses {
                continue;
            }
            winding_delta += if entry_y >= y_mid { -1 } else { 1 } * winding_factor;
        }
        if winding_delta != 0 {
            row_vec.push(BoundaryCellRecord { col: full_col, coverage: 0.0, winding_delta });
        }
    }
    Ok(())
}

/// Coverage-mode processing of one polygon.
pub(crate) fn process_polygon(
    poly: &Polygon,
    full_grid: &Grid<0>,
    poly_id: i32,
    runs: &mut Vec<GridRun>,
    edges: &mut Vec<GridEdge>,
) -> Result<(), String> {
    let Some(exterior) = poly.rings.first() else { return Ok(()) };
    if exterior.len() < 4 {
        return Ok(());
    }
    let dx = full_grid.dx();
    let dy = full_grid.dy();

    // Region of interest: exterior envelope clipped to the grid.
    let env = ring_envelope(exterior);
    if !env.intersects(full_grid.extent()) {
        return Ok(());
    }
    let region = full_grid.extent().intersection(&env);
    if region.is_empty() {
        return Ok(());
    }

    let subgrid_bounded = full_grid.shrink_to_fit(&region)?;
    let subgrid = subgrid_bounded.make_infinite();
    if subgrid.is_empty() {
        return Ok(());
    }
    let sub_rows = subgrid.rows() - 2;
    let sub_cols = subgrid.cols() - 2;
    let row_off = ((full_grid.ymax() - subgrid_bounded.ymax()) / dy).round() as usize;
    let col_off = ((subgrid_bounded.xmin() - full_grid.xmin()) / dx).round() as usize;

    let mut row_data: Vec<Vec<BoundaryCellRecord>> = vec![Vec::new(); sub_rows];
    let mut scratch = Scratch::default();

    for (r, ring) in poly.rings.iter().enumerate() {
        if ring.len() < 4 {
            continue;
        }
        walk_ring(
            ring.clone(),
            is_ccw(ring),
            r == 0,
            &subgrid,
            &mut row_data,
            sub_rows,
            sub_cols,
            col_off,
            &mut scratch,
        )?;
    }

    // Winding sweep per row.
    let tol: f32 = 1e-6;
    let one_minus_tol: f32 = 1.0 - tol;

    for (sr, row_vec) in row_data.iter_mut().enumerate() {
        if row_vec.is_empty() {
            continue;
        }
        row_vec.sort_by_key(|r| r.col);

        // Merge same-column entries in push order (stable sort).
        let mut merged: Vec<BoundaryCellRecord> = Vec::with_capacity(row_vec.len());
        for rec in row_vec.iter() {
            if let Some(last) = merged.last_mut() {
                if last.col == rec.col {
                    last.coverage += rec.coverage;
                    last.winding_delta += rec.winding_delta;
                    continue;
                }
            }
            merged.push(*rec);
        }

        let mut winding: i32 = 0;
        let mut prev_col: i64 = -2;
        let full_row = (row_off + sr) as i32;

        for mc in &merged {
            // Interior run between the previous boundary cell and this one.
            // prev_col > -2 means at least one cell (including padding) seen.
            if winding != 0 && prev_col > -2 && mc.col > prev_col + 1 {
                runs.push(GridRun {
                    row: full_row,
                    col_start: (prev_col + 1) as i32,
                    col_end: mc.col as i32,
                    id: poly_id,
                });
            }

            let w = mc.coverage;
            if w > tol && w < one_minus_tol {
                edges.push(GridEdge { row: full_row, col: mc.col as i32, fraction: w, id: poly_id });
            } else if w >= one_minus_tol {
                let c = mc.col as i32;
                runs.push(GridRun { row: full_row, col_start: c, col_end: c + 1, id: poly_id });
            }

            winding += mc.winding_delta;
            prev_col = mc.col;
        }
    }
    Ok(())
}
