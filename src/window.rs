//! Exact coverage and exact integrals for small rings in cell index space,
//! without the traversal tables.
//!
//! `burn` answers "which cells does this polygon cover, and by how much"
//! as sparse runs and edges, for polygons of any size. Some callers ask a
//! narrower question many times over: a resampler asks it once per
//! destination pixel (rwarp's exact footprint), a zonal mean asks only for
//! the weighted sum. This module answers those two cheaply.
//!
//! Both work in *cell index space*: column `i` spans `x` in `[i, i + 1]`,
//! row `j` spans `y` in `[j, j + 1]`, row 0 at the top (`y` grows down).
//! For a `GridSpec`, `x = (X - xmin) / dx` and `y = (ymax - Y) / dy`.
//!
//! Both rest on one identity. Split each ring edge at every cell wall, so
//! each piece lies in one cell `(i, j)`, with midpoint `xm` and height
//! `dy`. By Green's theorem, the area of the ring inside the part of row
//! `j` left of `x = k` is the sum of `(x - k) dy` over the pieces in that
//! row left of `k` (cell walls add nothing: on a vertical wall `x - k` is
//! zero or the wall is outside the region, on a horizontal wall `dy` is
//! zero). Differencing at `k = i` and `k = i + 1` gives
//!
//! ```text
//! coverage(i, j) = sum over pieces in (i, j) of (xm - (i + 1)) dy
//!                - sum over pieces in row j left of column i of dy
//! ```
//!
//! [`CoverageWindow`] books the first sum ("area") and the piece heights
//! ("cover") into a dense window and resolves the second with one prefix
//! pass per row. This is the cell accumulation of FreeType's "gray"
//! rasteriser and of font-rs, with the coverage kept exact rather than
//! quantised. There is no traversal chaining, no per-cell record and, once
//! the window is warm, no allocation.
//!
//! [`RowPrefix`] goes one step further for weighted sums. With `S(i, j)`
//! the sum of a field over the cells of row `j` left of column `i`, the
//! field's x-antiderivative inside cell `(i, j)` is
//! `F(x) = S(i, j) + (x - i) v(i, j)`, and the integral of the field over
//! a ring is the boundary integral of `F dy`. That costs one lookup per
//! piece: the work is proportional to the ring's perimeter in cells, not
//! its area, and an edge shared by two rings (neighbouring pixels of a
//! destination grid) is integrated once and added to one, subtracted from
//! the other. Summed-area tables (Crow 1984) do the same for axis-aligned
//! boxes; the boundary integral extends them to any ring.
//!
//! Both are exact for straight edges up to floating point, and agree with
//! `burn` in coverage mode to its `f32` fractions for simple rings. They
//! differ from `burn` in one respect: contributions are signed, so a
//! self-overlapping ring counts its overlap twice (or cancels it, for
//! opposite turning) instead of taking the non-zero winding rule.
//!
//! # Measures
//!
//! The identity holds for any area element `w(y) dx dy` that depends on the
//! row coordinate only: replace `(x - k) dy` with the integral of
//! `(x - k) w(y) dy` along the piece, which has a closed form because `x`
//! is linear in `y` along it. [`Measure::CosLatitude`] is the sphere for a
//! grid whose rows are parallels (any regular longitude/latitude grid):
//! `w(y) = cos(phi(y))`. [`CoverageWindow`] then reports each cell's
//! covered measure as a fraction of the cell's own, and
//! [`Measure::ring_area`] gives a ring's measure. The planar measure is the
//! default and its arithmetic is unchanged.
//!
//! C++ port note: both are header-only candidates for the C++
//! controlledburn (no dependencies, plain loops over `double`); the piece
//! walk below ports line for line.

/// How area is measured in cell index space; see the module docs.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum Measure {
    /// Area in index space: every cell has area 1.
    #[default]
    Planar,
    /// Area on the unit sphere per unit of `x`, for a grid whose rows are
    /// parallels: the area element is `cos(phi) dx dy` with latitude
    /// `phi = phi0 + dphi * y` in radians. For a regular longitude/latitude
    /// grid with column width `dlam` radians, multiply by `|dlam * dphi|`
    /// for steradians. `phi` must stay within the poles.
    CosLatitude { phi0: f64, dphi: f64 },
}

/// `sin(t) / t`.
#[inline]
fn sinc(t: f64) -> f64 {
    if t.abs() < 1e-4 {
        1.0 - t * t / 6.0
    } else {
        t.sin() / t
    }
}

/// `(sin z - z cos z) / z^3`, by its series for small `z`.
#[inline]
fn sin_minus_zcos_over_z3(z: f64) -> f64 {
    let z2 = z * z;
    if z.abs() < 0.5 {
        1.0 / 3.0 - z2 * (1.0 / 30.0 - z2 * (1.0 / 840.0 - z2 * (1.0 / 45360.0 - z2 / 3991680.0)))
    } else {
        (z.sin() - z * z.cos()) / (z2 * z)
    }
}

impl Measure {
    /// The sphere for rows of latitude `phi0 + dphi * y` (radians).
    pub fn cos_latitude(phi0: f64, dphi: f64) -> Self {
        Measure::CosLatitude { phi0, dphi }
    }

    /// For one straight piece with midpoint `(xm, ym)` and extent
    /// `(dx, dy)`: the integrals of `w(y) dy` ("cover") and of `x w(y) dy`
    /// ("area") along it.
    #[inline]
    fn piece(&self, xm: f64, ym: f64, dx: f64, dy: f64) -> (f64, f64) {
        match *self {
            Measure::Planar => (dy, xm * dy),
            Measure::CosLatitude { phi0, dphi } => {
                // With u = y - ym over [-h, h], h = dy / 2, and x = xm + s u:
                // int w du = 2 cos(phim) sin(c h) / c = cos(phim) dy sinc(c h),
                // int u w du = -sin(phim) 2 (sin z - z cos z) / c^2, z = c h.
                let phim = phi0 + dphi * ym;
                let h = 0.5 * dy;
                let z = dphi * h;
                let cover = phim.cos() * dy * sinc(z);
                let first = -phim.sin() * 2.0 * dphi * h * h * h * sin_minus_zcos_over_z3(z);
                (cover, xm * cover + (dx / dy) * first)
            }
        }
    }

    /// The measure of a cell in row `row` (`y` in `[row, row + 1]`): 1 for
    /// [`Measure::Planar`].
    pub fn cell_area(&self, row: i64) -> f64 {
        match *self {
            Measure::Planar => 1.0,
            Measure::CosLatitude { phi0, dphi } => (phi0 + dphi * (row as f64 + 0.5)).cos() * sinc(0.5 * dphi),
        }
    }

    /// Signed measure of a ring (closed implicitly), with the sign of
    /// [`ring_area`]. Exact for straight edges in index space.
    pub fn ring_area(&self, xs: &[f64], ys: &[f64]) -> f64 {
        if let Measure::Planar = self {
            return ring_area(xs, ys);
        }
        let n = xs.len().min(ys.len());
        let mut a = 0.0;
        for k in 0..n {
            let k1 = if k + 1 == n { 0 } else { k + 1 };
            let dy = ys[k1] - ys[k];
            if dy != 0.0 {
                let (xm, ym) = (0.5 * (xs[k] + xs[k1]), 0.5 * (ys[k] + ys[k1]));
                a += self.piece(xm, ym, xs[k1] - xs[k], dy).1;
            }
        }
        a
    }
}

/// Split segment `p -> q` at every cell wall and call `f(xm, ym, dx, dy)`
/// for each piece with non-zero height. Horizontal pieces carry neither
/// area nor cover and are skipped.
#[inline]
fn for_each_piece(p: (f64, f64), q: (f64, f64), mut f: impl FnMut(f64, f64, f64, f64)) {
    let ((x0, y0), (x1, y1)) = (p, q);
    let (dx, dy) = (x1 - x0, y1 - y0);
    if dy == 0.0 || !(dx.is_finite() && dy.is_finite()) {
        return;
    }
    // Next wall in each direction, as an integer-valued f64.
    let (mut kx, sx) = if dx > 0.0 {
        (x0.floor() + 1.0, 1.0)
    } else if dx < 0.0 {
        (x0.ceil() - 1.0, -1.0)
    } else {
        (0.0, 0.0)
    };
    let (mut ky, sy) = if dy > 0.0 { (y0.floor() + 1.0, 1.0) } else { (y0.ceil() - 1.0, -1.0) };
    let mut t = 0.0f64;
    loop {
        // Recomputed from the wall each step, so the parameters do not
        // drift along long segments.
        let tx = if sx != 0.0 { (kx - x0) / dx } else { f64::INFINITY };
        let ty = (ky - y0) / dy;
        let tn = tx.min(ty).min(1.0);
        if tn > t {
            let tm = 0.5 * (t + tn);
            f(x0 + tm * dx, y0 + tm * dy, (tn - t) * dx, (tn - t) * dy);
        }
        if tn >= 1.0 {
            break;
        }
        if tx <= tn {
            kx += sx;
        }
        if ty <= tn {
            ky += sy;
        }
        t = tn;
    }
}

/// Split segment `p -> q` at every cell wall inside the box `[x0, x1]` by
/// `[y0, y1]` (whole cells) and call `f(xm, ym, dx, dy)` for each piece, as
/// [`for_each_piece`] does, with the work bounded by the box rather than
/// the segment's length.
///
/// The parts of the segment above or below the box are dropped. A part
/// left or right of the box becomes a vertical segment on the box's left
/// or right wall, over the same rows: both callers book a piece by its
/// height (and, under a measure, its rows) and by its x clamped to the
/// box, so this changes no result.
/// Without it, a ring vertex thrown far away (a transform near a
/// singularity) costs time in proportion to its distance.
///
/// Each clip point lies on a box wall: that coordinate is the wall's,
/// exactly, and the other is interpolated from the nearer end, so a far
/// vertex costs no precision near the box.
#[inline]
fn for_each_piece_in(
    p: (f64, f64),
    q: (f64, f64),
    (x0, x1, y0, y1): (f64, f64, f64, f64),
    mut f: impl FnMut(f64, f64, f64, f64),
) {
    let (dx, dy) = (q.0 - p.0, q.1 - p.1);
    if dy == 0.0 || !(p.0.is_finite() && p.1.is_finite() && q.0.is_finite() && q.1.is_finite()) {
        return;
    }
    if p.1.max(q.1) <= y0 || p.1.min(q.1) >= y1 {
        return;
    }
    // The other coordinate where the segment's line meets x = x or y = y,
    // from the nearer end.
    let y_at = |x: f64| {
        let e = if (x - p.0).abs() <= (x - q.0).abs() { p } else { q };
        e.1 + (x - e.0) * (dy / dx)
    };
    let x_at = |y: f64| {
        let e = if (y - p.1).abs() <= (y - q.1).abs() { p } else { q };
        e.0 + (y - e.1) * (dx / dy)
    };
    // Rows: move each end that is above or below the box onto its wall.
    let clip_y = |e: (f64, f64)| {
        let y = e.1.clamp(y0, y1);
        if y == e.1 {
            e
        } else {
            (x_at(y), y)
        }
    };
    let (a, b) = (clip_y(p), clip_y(q));
    // Columns: split where a -> b crosses a side wall, in order along it, so
    // each part is wholly left of, inside or right of the box, then clamp
    // each part's ends to the box.
    let mut pts = [a, b, b, b];
    let mut n = 1;
    let walls = if b.0 >= a.0 { [x0, x1] } else { [x1, x0] };
    for xw in walls {
        if (a.0 - xw) * (b.0 - xw) < 0.0 {
            pts[n] = (xw, y_at(xw).clamp(a.1.min(b.1), a.1.max(b.1)));
            n += 1;
        }
    }
    pts[n] = b;
    for k in 0..n {
        let (u, v) = (pts[k], pts[k + 1]);
        for_each_piece((u.0.clamp(x0, x1), u.1), (v.0.clamp(x0, x1), v.1), &mut f);
    }
}

/// Signed area of a ring given as separate coordinate slices (closed
/// implicitly: the last point joins the first), `sum of x dy`. Positive
/// when the ring runs clockwise as drawn with `y` down (counter-clockwise
/// in a `y`-up frame).
pub fn ring_area(xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len().min(ys.len());
    let mut a = 0.0;
    for k in 0..n {
        let k1 = if k + 1 == n { 0 } else { k + 1 };
        a += 0.5 * (xs[k] + xs[k1]) * (ys[k1] - ys[k]);
    }
    a
}

/// Exact coverage of one or more rings over a dense window of cells in
/// index space. See the module docs.
///
/// ```
/// use controlledburn::CoverageWindow;
/// // A square from (0.5, 0.5) to (2.5, 2.5) covers a quarter of each
/// // corner cell, half of each edge cell and all of the centre.
/// let (xs, ys) = ([0.5, 2.5, 2.5, 0.5], [0.5, 0.5, 2.5, 2.5]);
/// let mut w = CoverageWindow::default();
/// assert!(w.reset_to_ring(&xs, &ys));
/// w.add_ring(&xs, &ys, 1.0);
/// let c = w.resolve();
/// assert!((c[0] - 0.25).abs() < 1e-12 && (c[1] - 0.5).abs() < 1e-12 && (c[4] - 1.0).abs() < 1e-12);
/// ```
#[derive(Debug, Default, Clone)]
pub struct CoverageWindow {
    col0: i64,
    row0: i64,
    ncol: usize,
    nrow: usize,
    area: Vec<f64>,
    cover: Vec<f64>,
    measure: Measure,
}

impl CoverageWindow {
    /// A window that measures area with `measure` (default
    /// [`Measure::Planar`]). Coverage is then each cell's covered measure
    /// as a fraction of the cell's.
    pub fn with_measure(measure: Measure) -> Self {
        CoverageWindow { measure, ..Default::default() }
    }

    /// Change the measure; it stays across [`reset`](Self::reset).
    pub fn set_measure(&mut self, measure: Measure) {
        self.measure = measure;
    }

    pub fn measure(&self) -> Measure {
        self.measure
    }

    /// Clear the window and place it over columns `col0..col0 + ncol` and
    /// rows `row0..row0 + nrow`. Reuses the buffers.
    pub fn reset(&mut self, col0: i64, row0: i64, ncol: usize, nrow: usize) {
        self.col0 = col0;
        self.row0 = row0;
        self.ncol = ncol;
        self.nrow = nrow;
        let n = ncol * nrow;
        self.area.clear();
        self.area.resize(n, 0.0);
        self.cover.clear();
        self.cover.resize(n, 0.0);
    }

    /// Clear the window and place it over every cell the ring's bounding
    /// box touches. Returns false (and leaves an empty window) when a
    /// coordinate is not finite.
    pub fn reset_to_ring(&mut self, xs: &[f64], ys: &[f64]) -> bool {
        let (mut x0, mut x1, mut y0, mut y1) = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
        for (&x, &y) in xs.iter().zip(ys) {
            x0 = x0.min(x);
            x1 = x1.max(x);
            y0 = y0.min(y);
            y1 = y1.max(y);
        }
        if !(x0.is_finite() && x1.is_finite() && y0.is_finite() && y1.is_finite()) {
            self.reset(0, 0, 0, 0);
            return false;
        }
        let (c0, r0) = (x0.floor() as i64, y0.floor() as i64);
        let c1 = (x1.ceil() as i64).max(c0 + 1);
        let r1 = (y1.ceil() as i64).max(r0 + 1);
        self.reset(c0, r0, (c1 - c0) as usize, (r1 - r0) as usize);
        true
    }

    pub fn col0(&self) -> i64 {
        self.col0
    }
    pub fn row0(&self) -> i64 {
        self.row0
    }
    pub fn ncol(&self) -> usize {
        self.ncol
    }
    pub fn nrow(&self) -> usize {
        self.nrow
    }

    /// Book one piece. A piece on the window's far wall (`x = col0 + ncol`
    /// exactly, the right edge of a bounding box) is booked to the last
    /// column; on a wall `x = k` a piece adds the same to cell `k` whether
    /// it is booked to cell `k - 1` or `k`, so this is exact as long as the
    /// area term uses the column it is booked to. Pieces outside the
    /// window's rows are dropped; pieces left or right of it clamp to its
    /// first or last column, which keeps the row sums right.
    ///
    /// Under a measure the area term is the integral of
    /// `(x - (i + 1)) w(y) dy` along the piece; a piece moved onto a wall by
    /// the clamp is vertical there.
    #[inline]
    fn book(&mut self, xm: f64, ym: f64, dx: f64, dy: f64, w: f64) {
        if self.ncol == 0 {
            return;
        }
        let lj = ym.floor() as i64 - self.row0;
        if lj < 0 || lj >= self.nrow as i64 {
            return;
        }
        let li = (xm.floor() as i64 - self.col0).clamp(0, self.ncol as i64 - 1);
        let i = (self.col0 + li) as f64;
        let x = xm.clamp(i, i + 1.0);
        let k = lj as usize * self.ncol + li as usize;
        if let Measure::Planar = self.measure {
            let dy = w * dy;
            self.area[k] += (x - (i + 1.0)) * dy;
            self.cover[k] += dy;
        } else {
            let dx = if x == xm { dx } else { 0.0 };
            let (cover, area) = self.measure.piece(x - (i + 1.0), ym, dx, dy);
            self.area[k] += w * area;
            self.cover[k] += w * cover;
        }
    }

    /// Add segment `p -> q` with weight `w` (orientation is the caller's).
    pub fn add_segment(&mut self, p: (f64, f64), q: (f64, f64), w: f64) {
        if self.ncol == 0 || self.nrow == 0 {
            return;
        }
        let (c0, r0) = (self.col0 as f64, self.row0 as f64);
        let bounds = (c0, c0 + self.ncol as f64, r0, r0 + self.nrow as f64);
        for_each_piece_in(p, q, bounds, |xm, ym, dx, dy| self.book(xm, ym, dx, dy, w));
    }

    /// Add a ring (closed implicitly) so that its interior counts `weight`
    /// whichever way it runs: 1 for an exterior, -1 for a hole.
    pub fn add_ring(&mut self, xs: &[f64], ys: &[f64], weight: f64) {
        let n = xs.len().min(ys.len());
        if n < 3 {
            return;
        }
        // The identity carries the sign of ring_area.
        let w = if ring_area(xs, ys) < 0.0 { -weight } else { weight };
        for k in 0..n {
            let k1 = if k + 1 == n { 0 } else { k + 1 };
            self.add_segment((xs[k], ys[k]), (xs[k1], ys[k1]), w);
        }
    }

    /// Resolve the booked pieces into coverage, row-major over the window,
    /// and return it. Call once, after the last `add_*`. Under a measure,
    /// coverage is the covered measure over the cell's.
    pub fn resolve(&mut self) -> &[f64] {
        for j in 0..self.nrow {
            let row = j * self.ncol..(j + 1) * self.ncol;
            let mut acc = 0.0;
            for (a, c) in self.area[row.clone()].iter_mut().zip(&self.cover[row.clone()]) {
                *a -= acc;
                acc += *c;
            }
            if !matches!(self.measure, Measure::Planar) {
                let cell = self.measure.cell_area(self.row0 + j as i64);
                for a in &mut self.area[row] {
                    *a /= cell;
                }
            }
        }
        &self.area
    }

    /// The resolved coverage (after [`resolve`](Self::resolve)).
    pub fn coverage(&self) -> &[f64] {
        &self.area
    }
}

/// Row prefix sums of `K` fields over an `ncol` by `nrow` grid of cells,
/// for exact integrals over rings by their boundary. See the module docs.
///
/// Typical use is `K = 2`: values (0 where missing) and a validity mask,
/// so that `integral[0] / integral[1]` is the area-weighted mean over the
/// valid part of a ring. Fields are zero outside the grid.
#[derive(Debug, Default, Clone)]
pub struct RowPrefix<const K: usize> {
    ncol: usize,
    nrow: usize,
    /// `nrow` rows of `ncol + 1` prefix sums, the first 0.
    sums: Vec<[f64; K]>,
}

impl<const K: usize> RowPrefix<K> {
    /// Build from `f(col, row)`, called in row-major order. Reuses the
    /// buffer.
    pub fn build(&mut self, ncol: usize, nrow: usize, mut f: impl FnMut(usize, usize) -> [f64; K]) {
        self.ncol = ncol;
        self.nrow = nrow;
        self.sums.clear();
        self.sums.reserve(nrow * (ncol + 1));
        for j in 0..nrow {
            let mut acc = [0.0; K];
            self.sums.push(acc);
            for i in 0..ncol {
                let v = f(i, j);
                for k in 0..K {
                    acc[k] += v[k];
                }
                self.sums.push(acc);
            }
        }
    }

    pub fn from_fn(ncol: usize, nrow: usize, f: impl FnMut(usize, usize) -> [f64; K]) -> Self {
        let mut p = RowPrefix { ncol: 0, nrow: 0, sums: Vec::new() };
        p.build(ncol, nrow, f);
        p
    }

    pub fn ncol(&self) -> usize {
        self.ncol
    }
    pub fn nrow(&self) -> usize {
        self.nrow
    }

    /// The grid as a box for [`for_each_piece_in`].
    fn bounds(&self) -> (f64, f64, f64, f64) {
        (0.0, self.ncol as f64, 0.0, self.nrow as f64)
    }

    /// `F dy` for one piece. Inside cell `(i, j)` the antiderivative is
    /// `S(i) + (x - i) v(i)` with `v(i) = S(i + 1) - S(i)`.
    #[inline]
    fn piece(&self, xm: f64, ym: f64, dy: f64, out: &mut [f64; K]) {
        let j = ym.floor();
        if j < 0.0 || j >= self.nrow as f64 {
            return;
        }
        let o = j as usize * (self.ncol + 1);
        let i = xm.floor();
        if i < 0.0 {
            return;
        }
        if i >= self.ncol as f64 {
            let s = &self.sums[o + self.ncol];
            for k in 0..K {
                out[k] += s[k] * dy;
            }
            return;
        }
        let iu = i as usize;
        let (s0, s1) = (&self.sums[o + iu], &self.sums[o + iu + 1]);
        let f = xm - i;
        for k in 0..K {
            out[k] += (s0[k] + f * (s1[k] - s0[k])) * dy;
        }
    }

    /// Boundary integral of `F dy` along an open polyline. For rings that
    /// share edges, integrate each shared edge once and add it to one ring
    /// and subtract it from the other.
    pub fn integrate_polyline(&self, xs: &[f64], ys: &[f64]) -> [f64; K] {
        let mut out = [0.0; K];
        for k in 1..xs.len().min(ys.len()) {
            for_each_piece_in((xs[k - 1], ys[k - 1]), (xs[k], ys[k]), self.bounds(), |xm, ym, _, dy| {
                self.piece(xm, ym, dy, &mut out)
            });
        }
        out
    }

    /// Integral of each field over a ring (closed implicitly), positive
    /// whichever way the ring runs.
    pub fn integrate_ring(&self, xs: &[f64], ys: &[f64]) -> [f64; K] {
        let n = xs.len().min(ys.len());
        let mut out = [0.0; K];
        if n < 3 {
            return out;
        }
        for k in 0..n {
            let k1 = if k + 1 == n { 0 } else { k + 1 };
            for_each_piece_in((xs[k], ys[k]), (xs[k1], ys[k1]), self.bounds(), |xm, ym, _, dy| {
                self.piece(xm, ym, dy, &mut out)
            });
        }
        // The boundary integral carries the sign of ring_area.
        if ring_area(xs, ys) < 0.0 {
            for v in out.iter_mut() {
                *v = -*v;
            }
        }
        out
    }
}
