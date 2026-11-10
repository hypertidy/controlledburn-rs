//! Analytical coverage fraction for single-traversal cells.
//!
//! For one edge traversal through a grid cell, the covered area (to the
//! left of the traversal for a CCW ring) is a simple polygon: the
//! traversal path plus the cell corners lying on the clockwise perimeter
//! arc from exit back to entry. This avoids the chain-chasing algorithm
//! in `ee::traversal_areas` for the common case.

use crate::ee::bbox::BBox;
use crate::ee::measures::{area_signed, area_signed_by};
use crate::ee::perimeter::perimeter_distance;
use crate::geometry::Coord;

/// Exact coverage fraction for a single traversal (entry at
/// `coords[0]`, exit at `coords[last]`, both on the cell boundary).
pub fn analytical_covered_fraction(bbox: &BBox, coords: &[Coord]) -> f64 {
    let cell_area = bbox.area();
    if cell_area <= 0.0 || coords.len() < 2 {
        return 0.0;
    }
    let perim = bbox.perimeter();

    let exit_pd = perimeter_distance(bbox, &coords[coords.len() - 1]);
    let entry_pd = perimeter_distance(bbox, &coords[0]);

    // CW distance from exit to entry (going backward along the perimeter).
    let arc = if exit_pd > entry_pd + 1e-12 {
        exit_pd - entry_pd
    } else if entry_pd > exit_pd + 1e-12 {
        perim - entry_pd + exit_pd
    } else {
        // Entry ~= exit (degenerate): the traversal starts and ends at the
        // same point; its own path is the polygon.
        let m = coords.len();
        let n = if coords[0] != coords[m - 1] { m + 1 } else { m };
        return area_signed_by(n, |i| if i < m { coords[i] } else { coords[0] }).abs() / cell_area;
    };

    let h = bbox.height();
    let w = bbox.width();
    let corner_coord = bbox.corners_ccw_from_bottom_left();
    let corner_pd = [0.0, h, h + w, 2.0 * h + w];

    let cw_from_exit = |pd: f64| -> f64 {
        let mut d = exit_pd - pd;
        if d < 0.0 {
            d += perim;
        }
        d
    };

    // Corners strictly inside the CW arc from exit to entry, tagged with
    // their CW distance from the exit.
    let mut in_arc: [(Coord, f64); 4] = [(Coord::default(), 0.0); 4];
    let mut n_in_arc = 0;
    for i in 0..4 {
        let d = cw_from_exit(corner_pd[i]);
        if d > 1e-12 && d < arc - 1e-12 {
            in_arc[n_in_arc] = (corner_coord[i], d);
            n_in_arc += 1;
        }
    }
    // Nearest first. std::sort on <= 4 elements; distances are distinct
    // for distinct corners so stability is irrelevant.
    in_arc[..n_in_arc].sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    // The polygon is the traversal, then the corners, then back to the
    // start; summed in place rather than collected (one allocation per
    // boundary cell otherwise).
    let m = coords.len();
    let at = |i: usize| {
        if i < m {
            coords[i]
        } else if i < m + n_in_arc {
            in_arc[i - m].0
        } else {
            coords[0]
        }
    };
    area_signed_by(m + n_in_arc + 1, at).abs() / cell_area
}

/// Coverage fraction for a closed ring lying entirely within one cell.
pub fn closed_ring_covered_fraction(bbox: &BBox, ring: &[Coord]) -> f64 {
    let cell_area = bbox.area();
    if cell_area <= 0.0 {
        return 0.0;
    }
    area_signed(ring).abs() / cell_area
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn half_cell() {
        let b = BBox::new(0., 0., 1., 1.);
        // Enter bottom at (0.5, 0), exit top at (0.5, 1); CCW ring means
        // interior is to the left (west): half the cell.
        let t = vec![Coord::new(0.5, 0.), Coord::new(0.5, 1.)];
        assert!((analytical_covered_fraction(&b, &t) - 0.5).abs() < 1e-12);
        // Diagonal corner cut: enter left at (0, 0.5), exit bottom at
        // (0.5, 0); the interior (left of travel) is everything but the
        // bottom-left corner triangle.
        let t = vec![Coord::new(0., 0.5), Coord::new(0.5, 0.)];
        assert!((analytical_covered_fraction(&b, &t) - 0.875).abs() < 1e-12);
        let t = vec![Coord::new(0.5, 0.), Coord::new(0., 0.5)];
        assert!((analytical_covered_fraction(&b, &t) - 0.125).abs() < 1e-12);
    }
}
