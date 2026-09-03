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
pub(crate) fn walk_polyline(mut coords: Vec<Coord>, grid: &Grid<1>, closed: bool) -> Result<Cells, String> {
    let mut cells = Cells::default();
    if coords.is_empty() {
        return Ok(cells);
    }
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
                arena.push(x.coord);
                trav.len += 1;
                trav.exit = Some(x.side);
                if x.coord != next {
                    last_exit = Some(x.coord);
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
