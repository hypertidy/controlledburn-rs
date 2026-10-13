//! The cell-stepping walker: follows a polyline through the padded grid,
//! recording per-cell traversals (entry side, coordinates, exit side).
//!
//! This is a port of controlledburn's `walk_polyline`, itself a rewrite of
//! exactextract's `RasterCellIntersection::process_line` without the
//! `Cell` class. The control flow is kept identical, including the
//! "incomplete initial traversal" re-append for closed rings, so that the
//! set of traversals per cell (and hence every coverage fraction) matches
//! the C++ bit for bit.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

use crate::ee::bbox::BBox;
use crate::ee::side::Side;
use crate::geometry::Coord;
use crate::grid::Grid;

/// One pass of the polyline through one cell. Coordinates live in the
/// walk's shared arena (`Cells::coords`) as the range `start..start + len`.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Traversal {
    start: usize,
    len: usize,
    pub entry: Option<Side>,
    pub exit: Option<Side>,
}

impl Traversal {
    #[inline]
    pub fn coords<'a>(&self, arena: &'a [Coord]) -> &'a [Coord] {
        &arena[self.start..self.start + self.len]
    }

    #[inline]
    pub fn traversed(&self) -> bool {
        self.entry.is_some() && self.exit.is_some()
    }

    #[inline]
    pub fn is_closed_ring(&self, arena: &[Coord]) -> bool {
        let c = self.coords(arena);
        c.len() >= 3 && c[0] == c[c.len() - 1]
    }

    pub fn multiple_unique_coordinates(&self, arena: &[Coord]) -> bool {
        let c = self.coords(arena);
        c[1..].iter().any(|x| *x != c[0])
    }
}

/// Per-cell traversal data for one cell (row, col) of the padded grid.
#[derive(Clone, Debug)]
pub(crate) struct CellRecord {
    pub row: usize,
    pub col: usize,
    pub bbox: BBox,
    pub traversals: Vec<Traversal>,
}

/// Cells visited by a walk, in first-visit order, plus the coordinate
/// arena all traversals index into. Consumers that need a deterministic
/// order sort by (row, col).
#[derive(Clone, Debug, Default)]
pub(crate) struct Cells {
    pub records: Vec<CellRecord>,
    pub coords: Vec<Coord>,
}

/// Minimal multiplicative hasher for the (row, col) index. The walker
/// mostly revisits the current cell or a neighbour, so the map stays hot.
#[derive(Default)]
struct CellHasher(u64);

impl Hasher for CellHasher {
    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        for b in bytes {
            self.write_u64(*b as u64);
        }
    }
    #[inline]
    fn write_u64(&mut self, v: u64) {
        self.0 = (self.0.rotate_left(5) ^ v).wrapping_mul(0x517c_c1b7_2722_0a95);
    }
    #[inline]
    fn write_usize(&mut self, v: usize) {
        self.write_u64(v as u64);
    }
}

type CellIndex = HashMap<(usize, usize), usize, BuildHasherDefault<CellHasher>>;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Location {
    Inside,
    Boundary,
    Outside,
}

#[inline]
fn point_location(bbox: &BBox, c: &Coord) -> Location {
    if bbox.strictly_contains(c) {
        Location::Inside
    } else if bbox.contains(c) {
        Location::Boundary
    } else {
        Location::Outside
    }
}

/// Walk `coords` through `grid`, returning the visited cells.
///
/// `closed` selects how an initial traversal that began strictly inside a
/// cell (no entry side) is handled once it exits:
///   - `true` (polygon ring): its coordinates are appended to the input so
///     the start cell is re-entered with a proper entry side on the second
///     pass, stitching the start back to the closing edge.
///   - `false` (open polyline): the partial traversal is kept as is; the
///     line consumer accepts `entry == None` as valid.
///
/// Errors correspond to the C++ `std::runtime_error` throws ("Never get
/// here") and are reported by the caller as a note.
/// Snap vertices that lie within a tiny tolerance of a cell edge exactly
/// onto that edge, using the same arithmetic `Grid::cell` uses to place
/// the edges.
///
/// A vertex parsed from input (0.3) and an edge computed from the grid
/// (3 * 0.1 = 0.30000000000000004) can differ by a few ulp. Left unsnapped,
/// such a vertex is classified OUTSIDE the cell it visually sits on the
/// corner of, `BBox::crossing` extrapolates past the segment end, and the
/// walk leaves the grid along the extended line (the C++ core runs until
/// `bad_alloc`; here the padded-grid bounds check turns it into an error
/// and an empty polygon). Where the walk does close, a degenerate corner
/// exit makes `left_hand_area` over-count the cell. Tolerance is relative
/// to the cell size; only degenerate near-edge vertices move, by less than
/// 1e-9 of a cell.
fn snap_to_lattice(coords: &mut [Coord], grid: &Grid<1>) {
    let (dx, dy) = (grid.dx(), grid.dy());
    if dx <= 0.0 || dy <= 0.0 || dx.is_nan() || dy.is_nan() {
        return;
    }
    let (x0, x1, y0, y1) = (grid.xmin(), grid.xmax(), grid.ymin(), grid.ymax());
    // padded grids carry one padding row/col on each side
    let ncol = (grid.cols() - 2) as f64;
    let nrow = (grid.rows() - 2) as f64;
    let (tol_x, tol_y) = (dx * 1e-9, dy * 1e-9);
    for c in coords.iter_mut() {
        let kx = ((c.x - x0) / dx).round();
        if kx >= 0.0 && kx <= ncol {
            let bx = if kx == ncol { x1 } else { x0 + kx * dx };
            if (c.x - bx).abs() <= tol_x {
                c.x = bx;
            }
        }
        let ky = ((y1 - c.y) / dy).round();
        if ky >= 0.0 && ky <= nrow {
            let by = if ky == nrow { y0 } else { y1 - ky * dy };
            if (c.y - by).abs() <= tol_y {
                c.y = by;
            }
        }
    }
}

pub(crate) fn walk_polyline(mut coords: Vec<Coord>, grid: &Grid<1>, closed: bool) -> Result<Cells, String> {
    let mut cells = Cells::default();
    if coords.is_empty() {
        return Ok(cells);
    }
    snap_to_lattice(&mut coords, grid);
    let n_input = coords.len();
    // Every input vertex produces at most one arena entry plus one exit
    // crossing per cell change; reserve generously once.
    cells.coords.reserve(coords.len() * 2 + 8);
    let arena = &mut cells.coords;
    let records = &mut cells.records;
    let mut index = CellIndex::default();
    // Cache of the last cell looked up: consecutive visits are usually
    // the same cell or a neighbour.
    let mut last: Option<((usize, usize), usize)> = None;
    let n_rows = grid.rows();
    let n_cols = grid.cols();

    let mut pos: usize = 0;
    let mut row = grid.get_row(coords[0].y);
    let mut col = grid.get_column(coords[0].x);
    let mut last_exit: Option<Coord> = None;

    while pos < coords.len() {
        if row >= n_rows || col >= n_cols {
            return Err("walker stepped outside the padded grid".into());
        }
        let ci = match last {
            Some((k, i)) if k == (row, col) => i,
            _ => {
                let i = *index.entry((row, col)).or_insert_with(|| {
                    records.push(CellRecord { row, col, bbox: grid.cell(row, col), traversals: Vec::new() });
                    records.len() - 1
                });
                last = Some(((row, col), i));
                i
            }
        };
        let bbox = records[ci].bbox;

        let mut trav = Traversal { start: arena.len(), len: 0, entry: None, exit: None };

        while pos < coords.len() {
            let next = match last_exit {
                Some(c) => c,
                None => coords[pos],
            };

            if trav.len == 0 {
                // First coordinate for this traversal: enter the cell.
                trav.entry = bbox.side(&next);
                arena.push(next);
                trav.len += 1;
                if last_exit.is_some() {
                    last_exit = None;
                } else {
                    pos += 1;
                }
                continue;
            }

            if point_location(&bbox, &next) != Location::Outside {
                arena.push(next);
                trav.len += 1;
                if last_exit.is_some() {
                    last_exit = None;
                } else {
                    pos += 1;
                }
            } else {
                // Outside: compute the exit crossing. Use the previous
                // ORIGINAL vertex for robustness (same as Cell::take).
                let from = if pos > 0 { coords[pos - 1] } else { arena[arena.len() - 1] };
                let x = bbox.crossing(&from, &next).ok_or_else(|| String::from("Never get here."))?;
                let (mut exit_c, mut exit_s) = (x.coord, x.side);
                // Rounding guard (see snap_to_lattice): a crossing that lies
                // beyond the segment end can only come from a target that is
                // outside the box by rounding error. Clamp the target onto
                // the box and exit there instead of extrapolating.
                {
                    let d_next = (next.x - from.x).powi(2) + (next.y - from.y).powi(2);
                    let d_exit = (exit_c.x - from.x).powi(2) + (exit_c.y - from.y).powi(2);
                    if !exit_c.x.is_finite() || !exit_c.y.is_finite() || d_exit > d_next * (1.0 + 1e-9) {
                        exit_c = Coord::new(next.x.clamp(bbox.xmin, bbox.xmax), next.y.clamp(bbox.ymin, bbox.ymax));
                        exit_s = if next.x < bbox.xmin {
                            Side::Left
                        } else if next.x > bbox.xmax {
                            Side::Right
                        } else if next.y < bbox.ymin {
                            Side::Bottom
                        } else {
                            Side::Top
                        };
                    }
                }
                arena.push(exit_c);
                trav.len += 1;
                trav.exit = Some(exit_s);
                if exit_c != next {
                    last_exit = Some(exit_c);
                }
                break;
            }
        }

        // Force exit if stuck on the boundary (Cell::force_exit).
        if trav.exit.is_none() && trav.len > 0 {
            let last = arena[arena.len() - 1];
            if point_location(&bbox, &last) == Location::Boundary {
                trav.exit = bbox.side(&last);
            }
        }

        let exited = trav.exit.is_some();
        let incomplete = exited && trav.entry.is_none();

        if incomplete && closed {
            // A closed ring is reprocessed at most once per starting cell;
            // if the coordinate list keeps growing the walk is not converging.
            if coords.len() > 4 * n_input + 64 {
                return Err("ring walk did not converge (vertex on cell boundary?)".into());
            }
            coords.extend_from_slice(trav.coords(arena));
        }

        let exit_side = trav.exit;
        records[ci].traversals.push(trav);

        if exited {
            match exit_side {
                Some(Side::Top) => row = row.checked_sub(1).ok_or("walker stepped above the padded grid")?,
                Some(Side::Bottom) => row += 1,
                Some(Side::Left) => col = col.checked_sub(1).ok_or("walker stepped left of the padded grid")?,
                Some(Side::Right) => col += 1,
                None => {}
            }
        }
    }

    Ok(cells)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::GridSpec;

    #[test]
    fn diagonal_line_touches_expected_cells() {
        let g = GridSpec::new(0., 0., 4., 4., 4, 4).bounded().make_infinite();
        let line = vec![Coord::new(0.5, 0.4), Coord::new(3.5, 3.4)];
        let cells = walk_polyline(line, &g, false).unwrap();
        // A 45-degree line offset from the corners clips 7 cells
        // (padded indices: row 4..1, col 1..4).
        let mut keys: Vec<_> = cells.records.iter().map(|c| (c.row, c.col)).collect();
        keys.sort();
        assert_eq!(keys, vec![(1, 4), (2, 3), (2, 4), (3, 2), (3, 3), (4, 1), (4, 2)]);
        // start cell has entry None
        assert!(cells.records[0].traversals[0].entry.is_none());
    }

    #[test]
    fn closed_ring_reappends_start() {
        let g = GridSpec::new(0., 0., 10., 10., 10, 10).bounded().make_infinite();
        let ring = vec![
            Coord::new(2.5, 2.5),
            Coord::new(6.5, 2.5),
            Coord::new(6.5, 6.5),
            Coord::new(2.5, 6.5),
            Coord::new(2.5, 2.5),
        ];
        let cells = walk_polyline(ring, &g, true).unwrap();
        // The start cell (row 8, col 3 in padded coords) gets two
        // traversals: the incomplete one and the stitched one.
        let start = cells.records.iter().find(|c| (c.row, c.col) == (8, 3)).unwrap();
        assert_eq!(start.traversals.len(), 2);
        assert!(start.traversals[0].entry.is_none());
        assert!(start.traversals[1].traversed());
        assert_eq!(cells.records.len(), 16);
    }
}
