// Derived from exactextract traversal_areas.cpp
// Copyright (c) 2018-2019 ISciences, LLC. Apache License 2.0.

use super::bbox::BBox;
use super::measures::area;
use super::perimeter::perimeter_distance;
#[cfg(test)]
use super::perimeter::perimeter_distance_ccw;
use crate::geometry::Coord;

/// A chain is either one of the supplied traversals (index into
/// `coord_lists`) or a cell corner (index into the corner array).
#[derive(Clone, Copy, Debug)]
struct Chain {
    start: f64,
    stop: f64,
    /// 0..n_lists = traversal, n_lists..n_lists+4 = corner
    which: usize,
    visited: bool,
}

/// Reusable buffers for `left_hand_area_with`, to avoid allocating per cell.
#[derive(Default, Debug)]
pub struct Scratch {
    chains: Vec<Chain>,
    /// Chain indices ordered by (start ascending, index ascending).
    order: Vec<usize>,
    coords: Vec<Coord>,
}

/// Reference search: the candidate with the smallest counter-clockwise
/// perimeter distance from `current`'s exit, first index winning ties.
/// This is exactly exactextract's loop; kept for the differential test.
#[cfg(test)]
fn next_chain_naive(chains: &[Chain], current: usize, kill: usize, perimeter: f64) -> usize {
    let mut min = usize::MAX;
    let mut min_distance = f64::MAX;
    let stop = chains[current].stop;
    for (i, candidate) in chains.iter().enumerate() {
        if candidate.visited && i != kill {
            continue;
        }
        let distance = perimeter_distance_ccw(stop, candidate.start, perimeter);
        if distance < min_distance {
            min_distance = distance;
            min = i;
        }
    }
    min
}

/// Same result as the naive search using `order` (chains sorted by start,
/// then index). `perimeter_distance_ccw(stop, start)` is `stop - start`
/// for `start <= stop` and wraps otherwise, so the nearest candidate is
/// the eligible chain with the largest `start <= stop`, or failing that
/// the largest `start` overall; among equal starts the lowest index wins,
/// which is the first of the run in `order`.
#[inline]
fn next_chain_sorted(chains: &[Chain], order: &[usize], current: usize, kill: usize) -> usize {
    let stop = chains[current].stop;
    let eligible = |i: usize| !chains[i].visited || i == kill;

    // Position of the first chain with start > stop.
    let hi = order.partition_point(|&i| chains[i].start <= stop);

    // Scan backwards through start <= stop, then wrap to the top.
    let scan = |range: std::ops::Range<usize>| -> Option<usize> {
        let mut k = range.end;
        while k > range.start {
            k -= 1;
            let i = order[k];
            if !eligible(i) {
                continue;
            }
            // Within a run of equal starts the earliest index is the
            // winner: step back to the first eligible chain of the run.
            let start = chains[i].start;
            let mut best = i;
            let mut j = k;
            while j > range.start && chains[order[j - 1]].start == start {
                j -= 1;
                if eligible(order[j]) {
                    best = order[j];
                }
            }
            return Some(best);
        }
        None
    };
    scan(0..hi).or_else(|| scan(hi..order.len())).unwrap_or(usize::MAX)
}

/// Area of the region to the left of the given traversals inside `bbox`,
/// obtained by chasing chains around the perimeter: from each traversal's
/// exit, the next chain (another traversal or a corner) is the nearest
/// one counter-clockwise along the perimeter. Each traversal must start
/// and end on the box boundary.
#[allow(dead_code)]
pub fn left_hand_area(bbox: &BBox, coord_lists: &[&[Coord]]) -> f64 {
    left_hand_area_with(bbox, coord_lists, &mut Scratch::default())
}

/// `left_hand_area` with caller-provided scratch space.
pub fn left_hand_area_with(bbox: &BBox, coord_lists: &[&[Coord]], scratch: &mut Scratch) -> f64 {
    let corners = bbox.corners_ccw_from_bottom_left();
    let height = bbox.height();
    let width = bbox.width();
    let n_lists = coord_lists.len();

    let chains = &mut scratch.chains;
    chains.clear();
    for (i, coords) in coord_lists.iter().enumerate() {
        let start = perimeter_distance(bbox, &coords[0]);
        let stop = perimeter_distance(bbox, &coords[coords.len() - 1]);
        chains.push(Chain { start, stop, which: i, visited: false });
    }
    let corner_pd = [0.0, height, height + width, 2.0 * height + width];
    for (i, pd) in corner_pd.iter().enumerate() {
        chains.push(Chain { start: *pd, stop: *pd, which: n_lists + i, visited: false });
    }

    let order = &mut scratch.order;
    order.clear();
    order.extend(0..chains.len());
    order.sort_by(|&a, &b| {
        chains[a].start.partial_cmp(&chains[b].start).unwrap_or(std::cmp::Ordering::Equal).then(a.cmp(&b))
    });

    let mut sum = 0.0;
    let coords = &mut scratch.coords;
    for first in 0..chains.len() {
        // Corner chains (single coordinate) never start a loop.
        if chains[first].visited || chains[first].which >= n_lists {
            continue;
        }
        coords.clear();
        let mut chain = first;
        loop {
            chains[chain].visited = true;
            let which = chains[chain].which;
            if which < n_lists {
                coords.extend_from_slice(coord_lists[which]);
            } else {
                coords.push(corners[which - n_lists]);
            }
            chain = next_chain_sorted(chains, order, chain, first);
            if chain == first {
                break;
            }
        }
        coords.push(coords[0]);
        sum += area(coords);
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sorted_search_matches_naive() {
        // Deterministic pseudo-random chains with plenty of ties.
        let mut state: u64 = 12345;
        let mut rnd = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % 16) as f64
        };
        let perimeter = 16.0;
        for _ in 0..2000 {
            let n = 2 + (rnd() as usize % 8);
            let mut chains: Vec<Chain> =
                (0..n).map(|i| Chain { start: rnd(), stop: rnd(), which: i, visited: rnd() < 6.0 }).collect();
            let mut order: Vec<usize> = (0..n).collect();
            order.sort_by(|&a, &b| chains[a].start.partial_cmp(&chains[b].start).unwrap().then(a.cmp(&b)));
            for current in 0..n {
                for kill in 0..n {
                    chains[kill].visited = true;
                    let a = next_chain_naive(&chains, current, kill, perimeter);
                    let b = next_chain_sorted(&chains, &order, current, kill);
                    assert_eq!(a, b, "chains={chains:?} current={current} kill={kill}");
                }
            }
        }
    }

    #[test]
    fn diagonal_half() {
        let b = BBox::new(0., 0., 1., 1.);
        // Enter on the left at (0, 0.5), exit on the bottom at (0.5, 0):
        // left of the travel direction is everything but the bottom-left
        // corner triangle.
        let t = vec![Coord::new(0., 0.5), Coord::new(0.5, 0.)];
        assert!((left_hand_area(&b, &[&t]) - 0.875).abs() < 1e-12);
        // A second parallel diagonal (top to right) cuts off the top-right
        // corner triangle, leaving the band between them.
        let t2 = vec![Coord::new(0.5, 1.0), Coord::new(1.0, 0.5)];
        assert!((left_hand_area(&b, &[&t, &t2]) - 0.75).abs() < 1e-12);
    }
}
