//! Regular grid specification and the exactextract-derived grid indexing.
//!
//! `GridSpec` is the public description of the target raster. `Grid<PAD>`
//! is the internal indexing structure ported from exactextract's
//! `Grid<extent_tag>`: `PAD = 0` is the bounded grid, `PAD = 1` is the
//! "infinite" grid with a one-cell padding ring on every side whose outer
//! cells extend to the finite limits of `f64`. Only the subset of the
//! original used by the burn engine is ported.
//!
//! Row 0 is the TOP row (the `ymax` edge).

use crate::ee::bbox::BBox;
use crate::error::BurnError;

/// Regular grid specification. Cell size is derived:
/// `dx = (xmax - xmin) / ncol`, `dy = (ymax - ymin) / nrow`.
/// Row 1 is the top row (the `ymax` edge), matching raster convention.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GridSpec {
    pub xmin: f64,
    pub ymin: f64,
    pub xmax: f64,
    pub ymax: f64,
    pub ncol: u32,
    pub nrow: u32,
}

impl GridSpec {
    pub const fn new(xmin: f64, ymin: f64, xmax: f64, ymax: f64, ncol: u32, nrow: u32) -> Self {
        GridSpec { xmin, ymin, xmax, ymax, ncol, nrow }
    }

    #[inline]
    pub fn dx(&self) -> f64 {
        (self.xmax - self.xmin) / self.ncol as f64
    }

    #[inline]
    pub fn dy(&self) -> f64 {
        (self.ymax - self.ymin) / self.nrow as f64
    }

    #[inline]
    pub fn extent(&self) -> BBox {
        BBox::new(self.xmin, self.ymin, self.xmax, self.ymax)
    }

    pub fn validate(&self) -> Result<(), BurnError> {
        if self.ncol == 0 || self.nrow == 0 {
            return Err(BurnError::InvalidGrid("ncol and nrow must be positive".into()));
        }
        if self.ncol > i32::MAX as u32 || self.nrow > i32::MAX as u32 {
            return Err(BurnError::InvalidGrid("ncol and nrow must fit in i32".into()));
        }
        if !(self.xmin.is_finite() && self.ymin.is_finite() && self.xmax.is_finite() && self.ymax.is_finite()) {
            return Err(BurnError::InvalidGrid("extent must be finite".into()));
        }
        if self.xmax <= self.xmin || self.ymax <= self.ymin {
            return Err(BurnError::InvalidGrid("invalid extent: xmax must be > xmin, ymax must be > ymin".into()));
        }
        Ok(())
    }

    pub(crate) fn bounded(&self) -> Grid<0> {
        Grid::new(self.extent(), self.dx(), self.dy())
    }
}

/// Port of `exactextract::Grid<extent_tag>`; `PAD` is the padding
/// (0 = bounded_extent, 1 = infinite_extent).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Grid<const PAD: usize> {
    extent: BBox,
    dx: f64,
    dy: f64,
    num_rows: usize,
    num_cols: usize,
}

impl<const PAD: usize> Grid<PAD> {
    pub fn new(extent: BBox, dx: f64, dy: f64) -> Self {
        let rows = if extent.ymax > extent.ymin { ((extent.ymax - extent.ymin) / dy).round() as usize } else { 0 };
        let cols = if extent.xmax > extent.xmin { ((extent.xmax - extent.xmin) / dx).round() as usize } else { 0 };
        Grid { extent, dx, dy, num_rows: 2 * PAD + rows, num_cols: 2 * PAD + cols }
    }

    /// Column index for `x`. For the padded grid, values outside the extent
    /// map to the padding columns; `x == xmax` maps to the last real
    /// column. For the bounded grid the caller must have range-checked `x`
    /// (the C++ throws `out_of_range`; here it is a debug assertion and
    /// the index is clamped).
    pub fn get_column(&self, x: f64) -> usize {
        if PAD > 0 {
            if x < self.extent.xmin {
                return 0;
            }
            if x > self.extent.xmax {
                return self.num_cols - 1;
            }
            if x == self.extent.xmax {
                return self.num_cols - 2;
            }
        } else {
            debug_assert!(x >= self.extent.xmin && x <= self.extent.xmax, "x out of range");
            if x == self.extent.xmax {
                return self.num_cols - 1;
            }
        }
        let raw = PAD + ((x - self.extent.xmin) / self.dx).floor() as usize;
        raw.min(self.get_column(self.extent.xmax))
    }

    /// Row index for `y`; see `get_column`.
    pub fn get_row(&self, y: f64) -> usize {
        if PAD > 0 {
            if y > self.extent.ymax {
                return 0;
            }
            if y < self.extent.ymin {
                return self.num_rows - 1;
            }
            if y == self.extent.ymin {
                return self.num_rows - 2;
            }
        } else {
            debug_assert!(y >= self.extent.ymin && y <= self.extent.ymax, "y out of range");
            if y == self.extent.ymin {
                return self.num_rows - 1;
            }
        }
        let raw = PAD + ((self.extent.ymax - y) / self.dy).floor() as usize;
        raw.min(self.get_row(self.extent.ymin))
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.num_rows <= 2 * PAD && self.num_cols <= 2 * PAD
    }

    #[inline]
    pub fn rows(&self) -> usize {
        self.num_rows
    }

    #[inline]
    pub fn cols(&self) -> usize {
        self.num_cols
    }

    #[inline]
    pub fn xmin(&self) -> f64 {
        self.extent.xmin
    }
    #[inline]
    pub fn xmax(&self) -> f64 {
        self.extent.xmax
    }
    #[inline]
    pub fn ymin(&self) -> f64 {
        self.extent.ymin
    }
    #[inline]
    pub fn ymax(&self) -> f64 {
        self.extent.ymax
    }
    #[inline]
    pub fn dx(&self) -> f64 {
        self.dx
    }
    #[inline]
    pub fn dy(&self) -> f64 {
        self.dy
    }
    #[inline]
    pub fn extent(&self) -> &BBox {
        &self.extent
    }

    /// Shrink the grid to the smallest sub-grid (snapped to this grid's
    /// cell boundaries) that contains `b`. Verbatim port of exactextract,
    /// including its floating-point corrections. `b` must lie within the
    /// extent.
    pub fn shrink_to_fit(&self, b: &BBox) -> Result<Self, String> {
        let e = &self.extent;
        if b.xmin < e.xmin || b.ymin < e.ymin || b.xmax > e.xmax || b.ymax > e.ymax {
            return Err("Cannot shrink extent to bounds larger than original.".into());
        }

        let mut col0 = self.get_column(b.xmin);
        let mut row1 = self.get_row(b.ymax);

        // Snap xmin and ymax to the upper-left corner of the supplied extent.
        let mut snapped_xmin = e.xmin + (col0 - PAD) as f64 * self.dx;
        let mut snapped_ymax = e.ymax - (row1 - PAD) as f64 * self.dy;

        // Because of floating-point round-off the snapped corner may not
        // actually contain the requested corner.
        if b.xmin < snapped_xmin {
            snapped_xmin -= self.dx;
            col0 -= 1;
        }
        if b.ymax > snapped_ymax {
            snapped_ymax += self.dy;
            row1 -= 1;
        }

        let col1 = self.get_column(b.xmax);
        let row0 = self.get_row(b.ymin);

        let mut num_rows = 1 + (row0 - row1);
        let mut num_cols = 1 + (col1 - col0);

        // If xmax or ymin falls cleanly on a cell boundary we need one
        // fewer row/column, because the rightmost cell is a closed interval
        // in x and the lowermost cell a closed interval in y.
        if num_rows > 2 && (snapped_ymax - (num_rows - 1) as f64 * self.dy <= b.ymin) {
            num_rows -= 1;
        }
        if num_cols > 2 && (snapped_xmin + (num_cols - 1) as f64 * self.dx >= b.xmax) {
            num_cols -= 1;
        }

        // Offsets relative to the new xmin/ymax origin, so that repeated
        // shrink calls with the same inputs give the same result.
        let mut reduced = BBox::new(
            snapped_xmin,
            (snapped_ymax - num_rows as f64 * self.dy).min(b.ymin),
            (snapped_xmin + num_cols as f64 * self.dx).max(b.xmax),
            snapped_ymax,
        );

        // Fudge the computed xmax and ymin, if needed, so the extent does
        // not grow during a shrink.
        if reduced.xmax > e.xmax {
            if ((reduced.xmax - reduced.xmin) / self.dx).round() == ((e.xmax - reduced.xmin) / self.dx).round() {
                reduced.xmax = e.xmax;
            } else {
                return Err("Shrink operation failed.".into());
            }
        }
        if reduced.ymin < e.ymin {
            if ((reduced.ymax - reduced.ymin) / self.dy).round() == ((reduced.ymax - e.ymin) / self.dy).round() {
                reduced.ymin = e.ymin;
            } else {
                return Err("Shrink operation failed.".into());
            }
        }

        let out = Grid::new(reduced, self.dx, self.dy);
        if b.xmin < out.xmin() || b.ymin < out.ymin() || b.xmax > out.xmax() || b.ymax > out.ymax() {
            return Err("Shrink operation failed.".into());
        }
        Ok(out)
    }
}

impl Grid<0> {
    /// Box of cell (row, col) in a bounded grid. Cells along the right and
    /// bottom edges are stretched to the extent so that floating-point
    /// error in `xmin + n*dx` never loses part of the extent.
    #[allow(dead_code)]
    pub fn cell(&self, row: usize, col: usize) -> BBox {
        BBox::new(
            self.xmin() + col as f64 * self.dx,
            if row == self.rows() - 1 { self.ymin() } else { self.ymax() - (row + 1) as f64 * self.dy },
            if col == self.cols() - 1 { self.xmax() } else { self.xmin() + (col + 1) as f64 * self.dx },
            self.ymax() - row as f64 * self.dy,
        )
    }

    pub fn make_infinite(&self) -> Grid<1> {
        Grid::new(self.extent, self.dx, self.dy)
    }
}

impl Grid<1> {
    /// Box of cell (row, col) in a padded grid. Padding cells extend to the
    /// finite limits of `f64`; the last real column/row's far edge is the
    /// extent's xmax/ymin rather than `xmin + n*dx`.
    pub fn cell(&self, row: usize, col: usize) -> BBox {
        let cols = self.cols();
        let rows = self.rows();

        let xmin = if col == 0 {
            f64::MIN
        } else if col == cols - 1 {
            self.xmax()
        } else {
            self.xmin() + (col - 1) as f64 * self.dx
        };
        let xmax = match cols - col {
            1 => f64::MAX,
            2 => self.xmax(),
            _ => self.xmin() + col as f64 * self.dx,
        };
        let ymax = if row == 0 {
            f64::MAX
        } else if row == rows - 1 {
            self.ymin()
        } else {
            self.ymax() - (row - 1) as f64 * self.dy
        };
        let ymin = match rows - row {
            1 => f64::MIN,
            2 => self.ymin(),
            _ => self.ymax() - row as f64 * self.dy,
        };
        BBox::new(xmin, ymin, xmax, ymax)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spec_validation() {
        assert!(GridSpec::new(0., 0., 10., 10., 10, 10).validate().is_ok());
        assert!(GridSpec::new(0., 0., 10., 10., 0, 10).validate().is_err());
        assert!(GridSpec::new(0., 0., 0., 10., 1, 10).validate().is_err());
        assert!(GridSpec::new(f64::NAN, 0., 10., 10., 1, 1).validate().is_err());
    }

    #[test]
    fn index_boundaries() {
        let g = GridSpec::new(0., 0., 10., 10., 10, 10).bounded();
        assert_eq!(g.get_column(0.0), 0);
        assert_eq!(g.get_column(10.0), 9);
        assert_eq!(g.get_column(9.9999), 9);
        assert_eq!(g.get_row(10.0), 0);
        assert_eq!(g.get_row(0.0), 9);
        assert_eq!(g.get_row(5.0), 5);

        let p = g.make_infinite();
        assert_eq!(p.cols(), 12);
        assert_eq!(p.get_column(-1.0), 0);
        assert_eq!(p.get_column(11.0), 11);
        assert_eq!(p.get_column(10.0), 10);
        assert_eq!(p.get_column(0.0), 1);
        assert_eq!(p.get_row(11.0), 0);
        assert_eq!(p.get_row(-1.0), 11);
        assert_eq!(p.get_row(0.0), 10);

        let c = p.cell(0, 0);
        assert_eq!((c.xmin, c.ymax), (f64::MIN, f64::MAX));
        let c = p.cell(1, 1);
        assert_eq!((c.xmin, c.ymin, c.xmax, c.ymax), (0.0, 9.0, 1.0, 10.0));
        let c = p.cell(10, 10);
        assert_eq!((c.xmin, c.ymin, c.xmax, c.ymax), (9.0, 0.0, 10.0, 1.0));
    }

    #[test]
    fn shrink() {
        let g = GridSpec::new(0., 0., 10., 10., 10, 10).bounded();
        let s = g.shrink_to_fit(&BBox::new(2.5, 4.5, 6.5, 8.5)).unwrap();
        assert_eq!((s.xmin(), s.ymin(), s.xmax(), s.ymax()), (2.0, 4.0, 7.0, 9.0));
        assert_eq!((s.rows(), s.cols()), (5, 5));
        // aligned: xmax on a boundary needs one fewer column
        let s = g.shrink_to_fit(&BBox::new(2.0, 4.0, 6.0, 8.0)).unwrap();
        assert_eq!((s.xmin(), s.ymin(), s.xmax(), s.ymax()), (2.0, 4.0, 6.0, 8.0));
        assert_eq!((s.rows(), s.cols()), (4, 4));
        assert!(g.shrink_to_fit(&BBox::new(-1.0, 0., 1., 1.)).is_err());
    }
}
