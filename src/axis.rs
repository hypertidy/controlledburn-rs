//! One axis of a grid: where its cell walls are, and the lookups the
//! engine needs from them.
//!
//! This is the seam for generalising the grid beyond a regular lattice
//! (rectilinear axes next; see the grid-geometry plan in the issue
//! tracker). Only the regular axis exists so far, and every method
//! reproduces, expression for expression, the arithmetic the engine used
//! before the seam was introduced, so output is bit-identical. Keep it
//! that way: the golden parity tests compare edge fractions as f32 bit
//! patterns.
//!
//! Positions along an axis are counted from its *origin* wall: `xmin` for
//! the x axis (columns ascend with x) and `ymax` for the y axis (row 0 is
//! the top row, so rows ascend as y descends).

/// Which way cell indices run along the axis.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Dir {
    /// Index 0 at `lo`; indices grow with the coordinate (the x axis).
    Ascending,
    /// Index 0 at `hi`; indices grow as the coordinate falls (the y axis).
    Descending,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Axis {
    /// `n` cells of width `d` between `lo` and `hi`.
    Regular { lo: f64, hi: f64, d: f64, n: usize, dir: Dir },
}

impl Axis {
    pub fn regular(lo: f64, hi: f64, d: f64, dir: Dir) -> Self {
        let n = if hi > lo { ((hi - lo) / d).round() as usize } else { 0 };
        Axis::Regular { lo, hi, d, n, dir }
    }

    /// Number of cells.
    #[inline]
    pub fn n(&self) -> usize {
        match *self {
            Axis::Regular { n, .. } => n,
        }
    }

    #[inline]
    pub fn lo(&self) -> f64 {
        match *self {
            Axis::Regular { lo, .. } => lo,
        }
    }

    #[inline]
    pub fn hi(&self) -> f64 {
        match *self {
            Axis::Regular { hi, .. } => hi,
        }
    }

    /// Cell width of a regular axis. Regular-only by nature; callers that
    /// use it are the places a rectilinear axis will need attention.
    #[inline]
    pub fn step(&self) -> f64 {
        match *self {
            Axis::Regular { d, .. } => d,
        }
    }

    /// Position of `v` in cell units, measured from the origin wall:
    /// `(v - lo) / d` ascending, `(hi - v) / d` descending. Its floor is
    /// the (unclamped) cell index.
    #[inline]
    pub fn units(&self, v: f64) -> f64 {
        match *self {
            Axis::Regular { lo, hi, d, dir, .. } => match dir {
                Dir::Ascending => (v - lo) / d,
                Dir::Descending => (hi - v) / d,
            },
        }
    }

    /// Wall `i` counted from the origin (wall 0 is the origin). Wall `n`
    /// is exactly the far bound, so `origin + n * d` round-off never loses
    /// part of the extent.
    #[inline]
    pub fn wall(&self, i: usize) -> f64 {
        match *self {
            Axis::Regular { lo, hi, d, n, dir } => match dir {
                Dir::Ascending => {
                    if i == n {
                        hi
                    } else {
                        lo + i as f64 * d
                    }
                }
                Dir::Descending => {
                    if i == n {
                        lo
                    } else {
                        hi - i as f64 * d
                    }
                }
            },
        }
    }

    /// Centre of cell `i`.
    #[inline]
    pub fn centre(&self, i: i64) -> f64 {
        match *self {
            Axis::Regular { lo, hi, d, dir, .. } => match dir {
                Dir::Ascending => lo + (i as f64 + 0.5) * d,
                Dir::Descending => hi - (i as f64 + 0.5) * d,
            },
        }
    }

    /// Move `v` exactly onto the nearest wall when it is within 1e-9 of a
    /// cell of it (see `snap_to_lattice` in the walker for why).
    #[inline]
    pub fn snap(&self, v: f64) -> f64 {
        match *self {
            Axis::Regular { d, n, .. } => {
                if d <= 0.0 || d.is_nan() {
                    return v;
                }
                let k = self.units(v).round();
                if k >= 0.0 && k <= n as f64 {
                    let b = self.wall(k as usize);
                    if (v - b).abs() <= d * 1e-9 {
                        return b;
                    }
                }
                v
            }
        }
    }

    /// The smallest run of whole cells covering `[b_lo, b_hi]` (which must
    /// lie within the axis), as a new axis plus the index of its first cell
    /// in this one. `first` is this axis's cell index containing the bound
    /// at the origin end (`b_lo` ascending, `b_hi` descending) and `last`
    /// the one containing the far-end bound, both before any adjustment.
    ///
    /// One axis's half of exactextract's `Grid::shrink_to_fit`, including
    /// its floating-point corrections; x and y never interact there, so
    /// splitting it per axis changes nothing.
    pub fn shrink(&self, first: usize, last: usize, b_lo: f64, b_hi: f64) -> Result<(Axis, usize), String> {
        match *self {
            Axis::Regular { lo, hi, d, dir, .. } => {
                let mut first = first;
                match dir {
                    Dir::Ascending => {
                        // Snap to the origin-side wall of the first cell.
                        let mut snapped = lo + first as f64 * d;
                        // Round-off may leave the snapped wall past the bound.
                        if b_lo < snapped {
                            snapped -= d;
                            first -= 1;
                        }
                        let mut num = 1 + (last - first);
                        // A bound exactly on a wall needs one fewer cell (the
                        // last cell is closed at its far wall).
                        if num > 2 && (snapped + (num - 1) as f64 * d >= b_hi) {
                            num -= 1;
                        }
                        let new_lo = snapped;
                        let mut new_hi = (snapped + num as f64 * d).max(b_hi);
                        // Don't let the extent grow during a shrink.
                        if new_hi > hi {
                            if ((new_hi - new_lo) / d).round() == ((hi - new_lo) / d).round() {
                                new_hi = hi;
                            } else {
                                return Err("Shrink operation failed.".into());
                            }
                        }
                        Ok((Axis::regular(new_lo, new_hi, d, dir), first))
                    }
                    Dir::Descending => {
                        let mut snapped = hi - first as f64 * d;
                        if b_hi > snapped {
                            snapped += d;
                            first -= 1;
                        }
                        let mut num = 1 + (last - first);
                        if num > 2 && (snapped - (num - 1) as f64 * d <= b_lo) {
                            num -= 1;
                        }
                        let new_hi = snapped;
                        let mut new_lo = (snapped - num as f64 * d).min(b_lo);
                        if new_lo < lo {
                            if ((new_hi - new_lo) / d).round() == ((new_hi - lo) / d).round() {
                                new_lo = lo;
                            } else {
                                return Err("Shrink operation failed.".into());
                            }
                        }
                        Ok((Axis::regular(new_lo, new_hi, d, dir), first))
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regular_lookups() {
        let x = Axis::regular(0.0, 10.0, 1.0, Dir::Ascending);
        let y = Axis::regular(0.0, 10.0, 1.0, Dir::Descending);
        assert_eq!((x.n(), y.n()), (10, 10));
        assert_eq!(x.units(2.5), 2.5);
        assert_eq!(y.units(2.5), 7.5);
        assert_eq!((x.wall(0), x.wall(3), x.wall(10)), (0.0, 3.0, 10.0));
        assert_eq!((y.wall(0), y.wall(3), y.wall(10)), (10.0, 7.0, 0.0));
        assert_eq!((x.centre(0), y.centre(0)), (0.5, 9.5));
    }

    #[test]
    fn snap_moves_only_near_wall_values() {
        let x = Axis::regular(0.0, 1.0, 0.1, Dir::Ascending);
        // 0.3 parsed vs wall 3 computed as 3 * 0.1
        assert_eq!(x.snap(0.3), 3.0 * 0.1);
        assert_eq!(x.snap(0.35), 0.35);
        assert_eq!(x.snap(1.0 - 1e-12), 1.0);
    }
}
