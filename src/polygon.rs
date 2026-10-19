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
    /// The ring (0 the exterior) the record comes from.
    ring: u32,
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
    ring: u32,
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
    let is_exterior = ring == 0;
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

        // A coverage record goes in for every in-grid cell the ring runs
        // through, even at zero coverage, so the sweep treats the cell as a
        // boundary cell rather than filling it as part of an interior run.
        // Edges lying along both walls of a cell (lattice-aligned input) can
        // leave it with zero coverage and zero net winding: one crossing on
        // each wall, cancelling. Without a record the sweep never saw the
        // cell and filled it with the winding from its left. A traversal
        // that only touches a corner (entry and exit a rounding error apart)
        // has no extent, and its zero coverage says nothing about the cell:
        // no record, the sweep fills it from the winding as before.
        let has_extent = frac == 0.0 && in_grid_cols && {
            let tol_x = (cr.bbox.xmax - cr.bbox.xmin) * 1e-9;
            let tol_y = (cr.bbox.ymax - cr.bbox.ymin) * 1e-9;
            valid.iter().any(|t| {
                let tc = t.coords(arena);
                tc.iter().skip(1).any(|p| (p.x - tc[0].x).abs() > tol_x || (p.y - tc[0].y).abs() > tol_y)
            })
        };
        if frac != 0.0 || has_extent {
            row_vec.push(BoundaryCellRecord {
                col: full_col,
                coverage: coverage_factor * frac,
                winding_delta: 0,
                ring,
            });
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
            // Half-open at the row centre: a point exactly on y_mid counts
            // as above it. A traversal that ends on the centre line and
            // the next one that starts there then register the crossing
            // exactly once, in one cell; with a strict test on both sides
            // neither did, and the row lost its winding from that column
            // on (crossings on a column wall at the row centre are common
            // once vertices sit on the lattice).
            let entry_above = entry_y >= y_mid;
            let exit_above = exit_y >= y_mid;
            if entry_above == exit_above {
                continue;
            }
            winding_delta += if entry_above { -1 } else { 1 } * winding_factor;
        }
        if winding_delta != 0 {
            row_vec.push(BoundaryCellRecord { col: full_col, coverage: 0.0, winding_delta, ring });
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
    let n_rings = poly.rings.len();
    let mut scratch = Scratch::default();

    for (r, ring) in poly.rings.iter().enumerate() {
        if ring.len() < 4 {
            continue;
        }
        walk_ring(
            ring.clone(),
            is_ccw(ring),
            r as u32,
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

    // Winding per ring, so a boundary cell gains the cover of rings that
    // enclose it without passing through it (the exterior around a hole's
    // boundary cell). Allocated once; each row resets only what it touched.
    let mut ring_winding: Vec<i32> = vec![0; n_rings];
    let mut touched: Vec<u32> = Vec::new();
    // Rings passing through each merged cell, flat, with a span per cell.
    let mut through: Vec<(u32, i32)> = Vec::new();
    let mut spans: Vec<(usize, usize)> = Vec::new();
    let sign = |r: u32| -> i32 {
        if r == 0 {
            1
        } else {
            -1
        }
    };

    for (sr, row_vec) in row_data.iter_mut().enumerate() {
        if row_vec.is_empty() {
            continue;
        }
        row_vec.sort_by_key(|r| r.col);

        // Merge same-column entries in push order (stable sort), keeping
        // the rings that pass through each cell (each ring once per cell).
        let mut merged: Vec<BoundaryCellRecord> = Vec::with_capacity(row_vec.len());
        through.clear();
        spans.clear();
        for rec in row_vec.iter() {
            if let Some(last) = merged.last_mut() {
                if last.col == rec.col {
                    last.coverage += rec.coverage;
                    last.winding_delta += rec.winding_delta;
                    let span = spans.last_mut().unwrap();
                    match through[span.0..span.1].iter_mut().find(|(r, _)| *r == rec.ring) {
                        Some((_, d)) => *d += rec.winding_delta,
                        None => {
                            through.push((rec.ring, rec.winding_delta));
                            span.1 += 1;
                        }
                    }
                    continue;
                }
            }
            merged.push(*rec);
            spans.push((through.len(), through.len() + 1));
            through.push((rec.ring, rec.winding_delta));
        }

        let mut winding: i32 = 0;
        // Sum of sign(r) over rings with non-zero winding: +1 inside the
        // exterior, -1 inside a hole. Kept up to date as windings change,
        // so a cell costs only the rings through it, not all the rings.
        let mut inside: i32 = 0;
        let mut prev_col: i64 = -2;
        let full_row = (row_off + sr) as i32;

        for (mc, &(s0, s1)) in merged.iter().zip(&spans) {
            let rings = &through[s0..s1];
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

            // +1 inside the exterior, -1 inside a hole, for rings that do
            // not pass through this cell
            let enclosing: i32 = inside
                - rings.iter().filter(|&&(t, _)| ring_winding[t as usize] != 0).map(|&(t, _)| sign(t)).sum::<i32>();
            let w = if enclosing != 0 { mc.coverage + enclosing as f32 } else { mc.coverage };
            if w > tol && w < one_minus_tol {
                edges.push(GridEdge { row: full_row, col: mc.col as i32, fraction: w, id: poly_id });
            } else if w >= one_minus_tol {
                let c = mc.col as i32;
                runs.push(GridRun { row: full_row, col_start: c, col_end: c + 1, id: poly_id });
            }

            winding += mc.winding_delta;
            prev_col = mc.col;
            for &(r, d) in rings {
                let before = ring_winding[r as usize];
                let after = before + d;
                if before == 0 && after != 0 {
                    inside += sign(r);
                    touched.push(r);
                } else if before != 0 && after == 0 {
                    inside -= sign(r);
                }
                ring_winding[r as usize] = after;
            }
        }
        for &r in &touched {
            ring_winding[r as usize] = 0;
        }
        touched.clear();
    }
    Ok(())
}
