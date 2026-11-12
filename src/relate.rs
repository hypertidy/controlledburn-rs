//! DE-9IM of polygons against grid cells, after the burn.
//!
//! A burn says, for every cell, how much of it a polygon covers. The
//! DE-9IM intersection matrix says how the polygon and the cell (a closed
//! box) meet: interiors, boundaries and exteriors, each against each. For
//! almost every cell the burn already answers that. A cell the polygon's
//! boundary does not reach is wholly inside (a run: [`IntersectionMatrix::COVERED`])
//! or wholly outside (no record: [`IntersectionMatrix::DISJOINT`]).
//! [`relate_cells`] computes the matrix for the remaining cells, those
//! whose closed box the boundary touches, from the geometry, and needs
//! only the segments near each of them.
//!
//! Method, per polygonal geometry:
//!
//! - **Candidates.** Each boundary segment is walked row by row; a cell is
//!   a candidate when the segment meets its closed box, decided with exact
//!   orientation tests. This includes cells the boundary only touches,
//!   along a wall or at a corner, which have no record.
//! - **BI, BB.** BI (boundary meets the open cell) is a strict separating
//!   axis test. BB is 1 when a segment runs along a wall, 0 when the
//!   boundary meets the outline only at points.
//! - **IB, EB.** The contacts split the cell's outline into open arcs, each
//!   wholly inside or outside the polygon. Going from one arc to the next
//!   past a contact, inside or outside flips once for each boundary branch
//!   that leaves the contact into the cell, which is a comparison of
//!   coordinates. One known arc per cell then fixes the rest: from a
//!   neighbour the boundary does not touch (its record is exact), or from
//!   a candidate neighbour already known through the shared wall. A group
//!   of candidates with neither (a polygon inside one cell, say) gets one
//!   exact point-in-polygon test.
//! - **II, EI.** Both 2 when the boundary crosses the open cell; otherwise
//!   the open cell is wholly inside or outside, as its arcs are.
//! - **IE, BE.** F only when the polygon lies inside the closed cell.
//!
//! The records alone would not do for candidates: the burn treats coverage
//! within 1e-6 of 0 or 1 as 0 or 1, so a cell holding a sliver of the
//! polygon, or of its exterior, has no record or a run. Off the boundary
//! that tolerance never applies, so the records are exact there.

use crate::approx::orient2d_sign;
use crate::burn::BurnOptions;
use crate::error::BurnError;
use crate::geometry::{Coord, Geometry, Polygon};
use crate::grid::GridSpec;
use crate::output::BurnResult;

/// A DE-9IM matrix: nine entries, row-major, the polygon's interior,
/// boundary and exterior (rows) against the cell's (columns). Each entry is
/// `F` (empty) or the dimension of the intersection, `0`, `1` or `2`.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct IntersectionMatrix([u8; 9]);

impl IntersectionMatrix {
    /// A cell wholly inside the polygon's interior (a run cell the
    /// boundary does not touch).
    pub const COVERED: IntersectionMatrix = IntersectionMatrix(*b"212FF1FF2");
    /// A cell the polygon does not touch.
    pub const DISJOINT: IntersectionMatrix = IntersectionMatrix(*b"FF2FF1212");

    /// The matrix as nine characters, e.g. `"212101212"`.
    pub fn as_str(&self) -> &str {
        // always ASCII
        std::str::from_utf8(&self.0).unwrap_or("?????????")
    }

    /// The entries as bytes, `b'F'`, `b'0'`, `b'1'` or `b'2'`.
    pub fn as_bytes(&self) -> &[u8; 9] {
        &self.0
    }

    /// Whether the matrix matches a DE-9IM pattern of nine characters from
    /// `T F * 0 1 2` (`T` is any of `0 1 2`, `*` anything). A pattern of
    /// another length never matches.
    pub fn matches(&self, pattern: &str) -> bool {
        let p = pattern.as_bytes();
        p.len() == 9
            && p.iter().zip(&self.0).all(|(&p, &v)| match p {
                b'*' => true,
                b'T' | b't' => v != b'F',
                b'F' | b'f' => v == b'F',
                _ => p == v,
            })
    }

    /// The polygon and the cell share a point.
    pub fn intersects(&self) -> bool {
        [0, 1, 3, 4].iter().any(|&i| self.0[i] != b'F')
    }

    /// The interiors meet: the polygon covers part of the cell's area.
    pub fn interiors_meet(&self) -> bool {
        self.0[0] != b'F'
    }

    /// The polygon covers the cell: no point of the cell is outside it.
    pub fn covers(&self) -> bool {
        self.0[6] == b'F' && self.0[7] == b'F'
    }

    /// They touch: they share a point, but not interior area.
    pub fn touches(&self) -> bool {
        self.0[0] == b'F' && self.intersects()
    }

    /// The cell covers the polygon: no point of the polygon is outside it.
    pub fn covered_by(&self) -> bool {
        self.0[2] == b'F' && self.0[5] == b'F'
    }
}

impl std::fmt::Display for IntersectionMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

impl std::fmt::Debug for IntersectionMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "IntersectionMatrix({})", self.as_str())
    }
}

/// The DE-9IM matrix of geometry `id` against cell `(row, col)`, for a
/// cell whose closed box the geometry's boundary touches.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CellRelation {
    pub row: i32,
    pub col: i32,
    pub id: i32,
    pub matrix: IntersectionMatrix,
}

/// DE-9IM matrices of polygons against the grid cells their boundaries
/// touch, from a finished burn and the same geometries.
///
/// `result` must be the burn of `geoms` onto `grid` with `opts` (pass the
/// same options; either mode, any `coalesce_runs` and `merge_parts`, and
/// `id_offset` gives the ids). The result may be cropped or sliced to a
/// window of the grid as long as the window's runs are all there; only
/// runs are read. Geometry `k` has `id = k + opts.id_offset`.
///
/// Returns one [`CellRelation`] for each in-grid cell whose closed box a
/// polygon's boundary meets, sorted by id, row and column. Every other
/// cell of a geometry is [`IntersectionMatrix::COVERED`] if it is in one of
/// that geometry's runs and [`IntersectionMatrix::DISJOINT`] if not, so
/// the two together give the matrix of every cell. Lines, points, empty
/// geometries and geometries with a non-finite coordinate give no
/// relations.
///
/// Cells are the burn's: column `c` spans `xmin + c * dx` to
/// `xmin + (c + 1) * dx`, the last column ending at `xmax` (rows the same,
/// from `ymax` down). Contacts are decided with exact predicates on the
/// input coordinates (no snapping); positions of crossings along a cell's
/// outline are computed in floating point, so two crossings of one wall
/// closer than its rounding count as one. Polygons are assumed valid
/// (OGC): rings that do not cross, holes inside their shell, parts of a
/// multi-polygon meeting only at points. On a non-affine grid burned in
/// index space, the matrix is the index-space one.
///
/// ```
/// # use controlledburn::{burn, relate_cells, BurnOptions, Coord, Geometry, GridSpec, IntersectionMatrix, Polygon};
/// let grid = GridSpec::new(0.0, 0.0, 4.0, 4.0, 4, 4);
/// // A square on cell walls, and a diamond whose corners touch them.
/// let square = Geometry::Polygon(Polygon::new(vec![vec![
///     Coord::new(1.0, 1.0), Coord::new(3.0, 1.0), Coord::new(3.0, 3.0), Coord::new(1.0, 3.0), Coord::new(1.0, 1.0),
/// ]]));
/// let geoms = [square];
/// let opts = BurnOptions::default();
/// let r = burn(&geoms, &grid, opts).unwrap();
/// let rel = relate_cells(&r, &geoms, &grid, opts).unwrap();
///
/// // The square's four cells are covered, and their boundary lies on walls.
/// let inside = rel.iter().filter(|c| c.matrix.covers()).count();
/// assert_eq!(inside, 4);
/// // The twelve cells around it touch it without sharing any area: the
/// // exact all-touched set is the 4 run cells plus these 12.
/// let touching = rel.iter().filter(|c| c.matrix.touches()).count();
/// assert_eq!(touching, 12);
/// assert!(rel.iter().all(|c| c.matrix != IntersectionMatrix::DISJOINT));
/// ```
pub fn relate_cells(
    result: &BurnResult,
    geoms: &[Geometry],
    grid: &GridSpec,
    opts: impl Into<BurnOptions>,
) -> Result<Vec<CellRelation>, BurnError> {
    let opts = opts.into();
    grid.validate()?;
    opts.validate()?;
    if !geoms.is_empty() {
        opts.id(geoms.len() - 1)?;
    }
    let runs = RunIndex::new(result);
    let walls = Walls::new(grid);
    let mut out = Vec::new();
    for (k, g) in geoms.iter().enumerate() {
        let id = opts.id(k)?;
        let polys: Vec<&Polygon> = match g {
            Geometry::Polygon(p) => vec![p],
            Geometry::MultiPolygon(ps) => ps.iter().collect(),
            _ => continue,
        };
        let Some(shape) = Shape::new(&polys) else { continue };
        relate_one(&shape, &walls, &runs, id, &mut out);
    }
    Ok(out)
}

// ------------------------------------------------------------------ cells

/// The burn's cell walls.
struct Walls {
    xmin: f64,
    xmax: f64,
    ymin: f64,
    ymax: f64,
    dx: f64,
    dy: f64,
    ncol: i64,
    nrow: i64,
}

impl Walls {
    fn new(g: &GridSpec) -> Self {
        Walls {
            xmin: g.xmin,
            xmax: g.xmax,
            ymin: g.ymin,
            ymax: g.ymax,
            dx: g.dx(),
            dy: g.dy(),
            ncol: g.ncol as i64,
            nrow: g.nrow as i64,
        }
    }

    /// As `Grid::cell`: the last column and row are stretched to the
    /// extent.
    fn cell(&self, row: i64, col: i64) -> CellBox {
        CellBox {
            x0: self.xmin + col as f64 * self.dx,
            x1: if col == self.ncol - 1 { self.xmax } else { self.xmin + (col + 1) as f64 * self.dx },
            y0: if row == self.nrow - 1 { self.ymin } else { self.ymax - (row + 1) as f64 * self.dy },
            y1: self.ymax - row as f64 * self.dy,
        }
    }

    /// Loose index of the column holding `x` (may be off by one at walls).
    fn col_of(&self, x: f64) -> i64 {
        ((x - self.xmin) / self.dx).floor() as i64
    }

    fn row_of(&self, y: f64) -> i64 {
        ((self.ymax - y) / self.dy).floor() as i64
    }
}

#[derive(Clone, Copy)]
struct CellBox {
    x0: f64,
    y0: f64,
    x1: f64,
    y1: f64,
}

/// The four walls in counterclockwise order from the bottom-left corner.
/// The outline is parameterised by `s`, its length from that corner.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Wall {
    Bottom,
    Right,
    Top,
    Left,
}

const WALLS: [Wall; 4] = [Wall::Bottom, Wall::Right, Wall::Top, Wall::Left];

impl CellBox {
    fn w(&self) -> f64 {
        self.x1 - self.x0
    }
    fn h(&self) -> f64 {
        self.y1 - self.y0
    }
    /// `s` at the start of each wall, and the perimeter.
    fn marks(&self) -> [f64; 5] {
        let (w, h) = (self.w(), self.h());
        let c2 = w + h;
        let c3 = c2 + w;
        [0.0, w, c2, c3, c3 + h]
    }

    /// Corners in counterclockwise order from the bottom left; corner `k`
    /// starts wall `k`.
    fn corners(&self) -> [Coord; 4] {
        [
            Coord::new(self.x0, self.y0),
            Coord::new(self.x1, self.y0),
            Coord::new(self.x1, self.y1),
            Coord::new(self.x0, self.y1),
        ]
    }

    fn strictly_inside(&self, p: Coord) -> bool {
        p.x > self.x0 && p.x < self.x1 && p.y > self.y0 && p.y < self.y1
    }

    /// `s` of a point on `wall`, given its coordinate along the wall.
    fn s_on(&self, wall: Wall, along: f64) -> f64 {
        let m = self.marks();
        match wall {
            Wall::Bottom => along - self.x0,
            Wall::Right => m[1] + (along - self.y0),
            Wall::Top => m[2] + (self.x1 - along),
            Wall::Left => m[3] + (self.y1 - along),
        }
    }

    /// The coordinate along `wall` of the point at `s`.
    fn along(&self, wall: Wall, s: f64) -> f64 {
        let m = self.marks();
        match wall {
            Wall::Bottom => self.x0 + s,
            Wall::Right => self.y0 + (s - m[1]),
            Wall::Top => self.x1 - (s - m[2]),
            Wall::Left => self.y1 - (s - m[3]),
        }
    }

    fn point_at(&self, s: f64) -> Coord {
        let m = self.marks();
        let s = s.rem_euclid(m[4]);
        if s < m[1] {
            Coord::new(self.along(Wall::Bottom, s), self.y0)
        } else if s < m[2] {
            Coord::new(self.x1, self.along(Wall::Right, s))
        } else if s < m[3] {
            Coord::new(self.along(Wall::Top, s), self.y1)
        } else {
            Coord::new(self.x0, self.along(Wall::Left, s))
        }
    }

    /// Whether a branch from a point on the open `wall` towards `o` goes
    /// into the open cell.
    fn inward(&self, wall: Wall, o: Coord) -> bool {
        match wall {
            Wall::Bottom => o.y > self.y0,
            Wall::Right => o.x < self.x1,
            Wall::Top => o.y < self.y1,
            Wall::Left => o.x > self.x0,
        }
    }

    /// Whether a branch from corner `k` towards `o` goes into the open cell.
    /// From a corner, the open cell is the open quadrant on its side.
    fn inward_corner(&self, k: usize, o: Coord) -> bool {
        let x_ok = if k == 0 || k == 3 { o.x > self.x0 } else { o.x < self.x1 };
        let y_ok = if k == 0 || k == 1 { o.y > self.y0 } else { o.y < self.y1 };
        x_ok && y_ok
    }
}

fn orient(a: Coord, b: Coord, c: Coord) -> i32 {
    orient2d_sign((a.x, a.y), (b.x, b.y), (c.x, c.y))
}

/// Segment `ab` meets the open box: strict separating axes, with the
/// box's corners on the segment's line counting as separated.
fn meets_open(a: Coord, b: Coord, bx: &CellBox) -> bool {
    if a.x.max(b.x) <= bx.x0 || a.x.min(b.x) >= bx.x1 || a.y.max(b.y) <= bx.y0 || a.y.min(b.y) >= bx.y1 {
        return false;
    }
    let o = bx.corners().map(|c| orient(a, b, c));
    !(o.iter().all(|&v| v >= 0) || o.iter().all(|&v| v <= 0))
}

/// Segment `ab` meets the closed box.
fn meets_closed(a: Coord, b: Coord, bx: &CellBox) -> bool {
    if a.x.max(b.x) < bx.x0 || a.x.min(b.x) > bx.x1 || a.y.max(b.y) < bx.y0 || a.y.min(b.y) > bx.y1 {
        return false;
    }
    let o = bx.corners().map(|c| orient(a, b, c));
    !(o.iter().all(|&v| v > 0) || o.iter().all(|&v| v < 0))
}

// ------------------------------------------------------------------ geometry

/// A polygonal geometry as one set of boundary segments.
struct Shape {
    segs: Vec<(Coord, Coord)>,
    bbox: [f64; 4],
}

impl Shape {
    fn new(polys: &[&Polygon]) -> Option<Shape> {
        let mut segs = Vec::new();
        let mut bbox = [f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY];
        for p in polys {
            for ring in &p.rings {
                if ring.len() < 3 {
                    continue;
                }
                for c in ring {
                    if !c.x.is_finite() || !c.y.is_finite() {
                        return None;
                    }
                    bbox = [bbox[0].min(c.x), bbox[1].min(c.y), bbox[2].max(c.x), bbox[3].max(c.y)];
                }
                let closed = ring.first() == ring.last();
                let n = ring.len();
                for i in 0..n - usize::from(closed) {
                    let (a, b) = (ring[i], ring[(i + 1) % n]);
                    if a != b {
                        segs.push((a, b));
                    }
                }
            }
        }
        (!segs.is_empty()).then_some(Shape { segs, bbox })
    }

    /// Exact point location by crossing number: true inside. `p` must not
    /// lie on the boundary.
    fn contains(&self, p: Coord) -> bool {
        let mut inside = false;
        for &(a, b) in &self.segs {
            if (a.y > p.y) != (b.y > p.y) {
                let o = orient(a, b, p);
                if (b.y > a.y && o > 0) || (b.y < a.y && o < 0) {
                    inside = !inside;
                }
            }
        }
        inside
    }
}

/// Runs by (id, row), merged, for "is this cell in a run of this id".
struct RunIndex {
    runs: Vec<(i32, i32, i32, i32)>,
}

impl RunIndex {
    fn new(r: &BurnResult) -> Self {
        let mut v: Vec<(i32, i32, i32, i32)> = r.runs.iter().map(|x| (x.id, x.row, x.col_start, x.col_end)).collect();
        v.sort_unstable();
        // merge overlapping or touching runs of one (id, row)
        let mut runs: Vec<(i32, i32, i32, i32)> = Vec::with_capacity(v.len());
        for x in v {
            match runs.last_mut() {
                Some(l) if l.0 == x.0 && l.1 == x.1 && x.2 <= l.3 => l.3 = l.3.max(x.3),
                _ => runs.push(x),
            }
        }
        RunIndex { runs }
    }

    fn covered(&self, id: i32, row: i32, col: i32) -> bool {
        let i = self.runs.partition_point(|&(a, r, c0, _)| (a, r, c0) <= (id, row, col));
        i > 0 && {
            let (a, r, _, c1) = self.runs[i - 1];
            a == id && r == row && col < c1
        }
    }
}

// ------------------------------------------------------------------ one geometry

/// A contact cluster event on the outline: a span `lo < hi` along a wall,
/// or a point `lo == hi`, with the number of branches leaving it into the
/// open cell.
#[derive(Clone, Copy)]
struct Event {
    lo: f64,
    hi: f64,
    inward: u32,
}

/// One candidate cell's outline: its open arcs as `(lo, hi)` in `s`, with
/// `hi` past the perimeter for an arc across the bottom-left corner, and
/// each arc's inside/outside relative to arc 0.
struct Outline {
    bx: CellBox,
    arcs: Vec<(f64, f64)>,
    parity: Vec<bool>,
    base: Option<bool>,
}

impl Outline {
    fn status(&self, j: usize) -> Option<bool> {
        self.base.map(|b| b ^ self.parity[j])
    }

    fn arc_at(&self, s: f64) -> Option<usize> {
        let per = self.bx.marks()[4];
        self.arcs.iter().position(|&(lo, hi)| (lo < s && s < hi) || (lo < s + per && s + per < hi))
    }

    /// A point strictly inside an arc on `wall`, as the coordinate along
    /// the wall, with the arc's index.
    fn on_wall(&self, wall: Wall) -> Option<(f64, usize)> {
        let m = self.bx.marks();
        let k = wall as usize;
        let (ws, we) = (m[k], m[k + 1]);
        for (j, &(lo, hi)) in self.arcs.iter().enumerate() {
            for shift in [0.0, m[4]] {
                let (a, b) = (lo.max(ws + shift) - shift, hi.min(we + shift) - shift);
                if a < b {
                    let along = self.bx.along(wall, (a + b) / 2.0);
                    return Some((along, j));
                }
            }
        }
        None
    }
}

struct CellFacts {
    bi: bool,
    bb: u8,
    outline: Outline,
}

fn relate_one(shape: &Shape, walls: &Walls, runs: &RunIndex, id: i32, out: &mut Vec<CellRelation>) {
    // candidate (row, col, segment) triples
    let mut hits: Vec<(i32, i32, u32)> = Vec::new();
    for (si, &(a, b)) in shape.segs.iter().enumerate() {
        let (miny, maxy) = (a.y.min(b.y), a.y.max(b.y));
        let r_lo = (walls.row_of(maxy) - 1).max(0);
        let r_hi = (walls.row_of(miny) + 1).min(walls.nrow - 1);
        for row in r_lo..=r_hi {
            let probe = walls.cell(row, 0);
            if maxy < probe.y0 || miny > probe.y1 {
                continue;
            }
            let (xa, xb) = if a.y == b.y {
                (a.x.min(b.x), a.x.max(b.x))
            } else {
                let x_at = |y: f64| a.x + (y - a.y) * (b.x - a.x) / (b.y - a.y);
                let (p, q) = (x_at(probe.y0.max(miny)), x_at(probe.y1.min(maxy)));
                (p.min(q), p.max(q))
            };
            let c_lo = (walls.col_of(xa) - 1).max(0);
            let c_hi = (walls.col_of(xb) + 1).min(walls.ncol - 1);
            for col in c_lo..=c_hi {
                if meets_closed(a, b, &walls.cell(row, col)) {
                    hits.push((row as i32, col as i32, si as u32));
                }
            }
        }
    }
    if hits.is_empty() {
        return;
    }
    hits.sort_unstable();
    hits.dedup();

    let mut cells: Vec<(i32, i32)> = Vec::new();
    let mut facts: Vec<CellFacts> = Vec::new();
    let mut i = 0;
    while i < hits.len() {
        let (row, col) = (hits[i].0, hits[i].1);
        let mut j = i;
        while j < hits.len() && (hits[j].0, hits[j].1) == (row, col) {
            j += 1;
        }
        let bx = walls.cell(row as i64, col as i64);
        facts.push(cell_facts(shape, &bx, hits[i..j].iter().map(|h| h.2 as usize)));
        cells.push((row, col));
        i = j;
    }

    // Inside/outside of each cell's arcs: seed from neighbours the boundary
    // does not touch, then across shared walls.
    let find = |r: i32, c: i32| cells.binary_search(&(r, c)).ok();
    let in_grid = |r: i32, c: i32| r >= 0 && c >= 0 && (r as i64) < walls.nrow && (c as i64) < walls.ncol;
    // neighbour across each wall, and the wall it sees
    let across = |w: Wall| -> (i32, i32, Wall) {
        match w {
            Wall::Bottom => (1, 0, Wall::Top),
            Wall::Right => (0, 1, Wall::Left),
            Wall::Top => (-1, 0, Wall::Bottom),
            Wall::Left => (0, -1, Wall::Right),
        }
    };
    let mut queue: Vec<usize> = Vec::new();
    for (k, &(row, col)) in cells.iter().enumerate() {
        let f = &mut facts[k];
        if f.outline.arcs.is_empty() {
            continue;
        }
        for w in WALLS {
            let (dr, dc, _) = across(w);
            let (r2, c2) = (row + dr, col + dc);
            if !in_grid(r2, c2) || find(r2, c2).is_some() {
                continue;
            }
            // No boundary on this wall: one arc covers its middle.
            let m = f.outline.bx.marks();
            let wi = w as usize;
            if let Some(j) = f.outline.arc_at((m[wi] + m[wi + 1]) / 2.0) {
                f.outline.base = Some(runs.covered(id, r2, c2) ^ f.outline.parity[j]);
                queue.push(k);
                break;
            }
        }
    }
    let mut next = 0;
    loop {
        while next < queue.len() {
            let k = queue[next];
            next += 1;
            let (row, col) = cells[k];
            for w in WALLS {
                let (dr, dc, w2) = across(w);
                let Some(k2) = find(row + dr, col + dc) else { continue };
                if facts[k2].outline.base.is_some() || facts[k2].outline.arcs.is_empty() {
                    continue;
                }
                let Some((along, j)) = facts[k].outline.on_wall(w) else { continue };
                let st = facts[k].outline.status(j).unwrap_or(false);
                let o2 = &mut facts[k2].outline;
                if let Some(j2) = o2.arc_at(o2.bx.s_on(w2, along)) {
                    o2.base = Some(st ^ o2.parity[j2]);
                    queue.push(k2);
                }
            }
        }
        // A group no neighbour reaches: one exact point-in-polygon test.
        let Some(k) = facts.iter().position(|f| f.outline.base.is_none() && !f.outline.arcs.is_empty()) else { break };
        let o = &mut facts[k].outline;
        let (lo, hi) = o.arcs[0];
        o.base = Some(shape.contains(o.bx.point_at((lo + hi) / 2.0)));
        queue.push(k);
    }

    for (k, &(row, col)) in cells.iter().enumerate() {
        let f = &facts[k];
        let o = &f.outline;
        let (mut ib, mut eb) = (false, false);
        for j in 0..o.arcs.len() {
            if o.status(j).unwrap_or(false) {
                ib = true;
            } else {
                eb = true;
            }
        }
        let (ii, ei) = if f.bi {
            (true, true)
        } else if o.arcs.is_empty() {
            // The whole outline is boundary (the polygon or a hole is this
            // cell): its coverage is exactly 0 or 1, so the record is exact.
            let c = runs.covered(id, row, col);
            (c, !c)
        } else {
            (ib, eb)
        };
        let bx = &o.bx;
        let within =
            shape.bbox[0] >= bx.x0 && shape.bbox[2] <= bx.x1 && shape.bbox[1] >= bx.y0 && shape.bbox[3] <= bx.y1;
        let d = |on: bool, dim: u8| if on { dim } else { b'F' };
        let m = [
            d(ii, b'2'),
            d(ib, b'1'),
            d(!within, b'2'),
            d(f.bi, b'1'),
            f.bb,
            d(!within, b'1'),
            d(ei, b'2'),
            d(eb, b'1'),
            b'2',
        ];
        out.push(CellRelation { row, col, id, matrix: IntersectionMatrix(m) });
    }
}

/// BI, BB and the outline's arcs for one candidate cell, from the
/// segments meeting its closed box.
fn cell_facts(shape: &Shape, bx: &CellBox, segs: impl Iterator<Item = usize>) -> CellFacts {
    let mut bi = false;
    let (mut line, mut point) = (false, false);
    let mut ev: Vec<Event> = Vec::new();
    let corners = bx.corners();
    let m = bx.marks();
    for si in segs {
        let (a, b) = shape.segs[si];
        bi |= meets_open(a, b, bx);
        if bx.strictly_inside(a) && bx.strictly_inside(b) {
            continue;
        }
        // Corners: the segment passes through or ends at one.
        for (k, &c) in corners.iter().enumerate() {
            if c.x < a.x.min(b.x) || c.x > a.x.max(b.x) || c.y < a.y.min(b.y) || c.y > a.y.max(b.y) {
                continue;
            }
            if orient(a, b, c) != 0 {
                continue;
            }
            point = true;
            let branches: &[Coord] = if c == a {
                &[b]
            } else if c == b {
                &[a]
            } else {
                &[a, b]
            };
            let inward = branches.iter().filter(|&&o| bx.inward_corner(k, o)).count() as u32;
            ev.push(Event { lo: m[k], hi: m[k], inward });
        }
        // Along a wall: a span, no branch into the cell.
        let mut along_wall = false;
        if a.y == b.y && (a.y == bx.y0 || a.y == bx.y1) {
            let (lo, hi) = (a.x.min(b.x).max(bx.x0), a.x.max(b.x).min(bx.x1));
            if lo < hi {
                let w = if a.y == bx.y0 { Wall::Bottom } else { Wall::Top };
                let (s0, s1) = (bx.s_on(w, lo), bx.s_on(w, hi));
                ev.push(Event { lo: s0.min(s1), hi: s0.max(s1), inward: 0 });
                along_wall = true;
            }
        }
        if a.x == b.x && (a.x == bx.x0 || a.x == bx.x1) {
            let (lo, hi) = (a.y.min(b.y).max(bx.y0), a.y.max(b.y).min(bx.y1));
            if lo < hi {
                let w = if a.x == bx.x1 { Wall::Right } else { Wall::Left };
                let (s0, s1) = (bx.s_on(w, lo), bx.s_on(w, hi));
                ev.push(Event { lo: s0.min(s1), hi: s0.max(s1), inward: 0 });
                along_wall = true;
            }
        }
        if along_wall {
            line = true;
            continue;
        }
        // Open walls: an endpoint on one, or a crossing between its corners.
        for w in WALLS {
            let (vertical, at, lo, hi, c_lo, c_hi) = match w {
                Wall::Bottom => (false, bx.y0, bx.x0, bx.x1, corners[0], corners[1]),
                Wall::Right => (true, bx.x1, bx.y0, bx.y1, corners[1], corners[2]),
                Wall::Top => (false, bx.y1, bx.x0, bx.x1, corners[3], corners[2]),
                Wall::Left => (true, bx.x0, bx.y0, bx.y1, corners[0], corners[3]),
            };
            let across = |p: Coord| if vertical { p.x } else { p.y };
            let along = |p: Coord| if vertical { p.y } else { p.x };
            for (e, o) in [(a, b), (b, a)] {
                if across(e) == at && along(e) > lo && along(e) < hi {
                    point = true;
                    ev.push(Event {
                        lo: bx.s_on(w, along(e)),
                        hi: bx.s_on(w, along(e)),
                        inward: u32::from(bx.inward(w, o)),
                    });
                }
            }
            let (da, db) = (across(a) - at, across(b) - at);
            if (da < 0.0 && db > 0.0) || (da > 0.0 && db < 0.0) {
                let (o1, o2) = (orient(a, b, c_lo), orient(a, b, c_hi));
                if o1 != 0 && o2 != 0 && o1 != o2 {
                    point = true;
                    let t = (at - across(a)) / (across(b) - across(a));
                    let pos = (along(a) + t * (along(b) - along(a))).clamp(lo, hi);
                    let s = bx.s_on(w, pos);
                    ev.push(Event { lo: s, hi: s, inward: 1 });
                }
            }
        }
    }
    let bb = if line {
        b'1'
    } else if point {
        b'0'
    } else {
        b'F'
    };
    let (arcs, parity) = arcs_of(ev, m[4]);
    CellFacts { bi, bb, outline: Outline { bx: *bx, arcs, parity, base: None } }
}

/// The open arcs between contact clusters, and each arc's inside/outside
/// relative to the first: crossing a cluster flips it once per branch the
/// cluster sends into the cell.
fn arcs_of(mut ev: Vec<Event>, per: f64) -> (Vec<(f64, f64)>, Vec<bool>) {
    if ev.is_empty() {
        return (vec![(0.0, per)], vec![false]);
    }
    ev.sort_by(|a, b| a.lo.total_cmp(&b.lo));
    let mut cl: Vec<Event> = Vec::new();
    for e in ev {
        match cl.last_mut() {
            Some(c) if e.lo <= c.hi => {
                c.hi = c.hi.max(e.hi);
                c.inward += e.inward;
            }
            _ => cl.push(e),
        }
    }
    // A cluster touching the end of the outline joins the one at its start.
    if cl.len() > 1 && cl[cl.len() - 1].hi >= per && cl[0].lo <= 0.0 {
        let last = cl.pop().unwrap_or(Event { lo: 0.0, hi: 0.0, inward: 0 });
        cl[0].lo = last.lo - per;
        cl[0].inward += last.inward;
    }
    if cl.len() == 1 && cl[0].hi - cl[0].lo >= per {
        return (Vec::new(), Vec::new());
    }
    let n = cl.len();
    let mut arcs = Vec::with_capacity(n);
    let mut parity = Vec::with_capacity(n);
    let mut p = false;
    for i in 0..n {
        let (lo, hi) = if i + 1 < n { (cl[i].hi, cl[i + 1].lo) } else { (cl[i].hi, cl[0].lo + per) };
        if i > 0 {
            p ^= cl[i].inward % 2 == 1;
        }
        if lo < hi {
            arcs.push((lo, hi));
            parity.push(p);
        }
    }
    (arcs, parity)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matrix_predicates() {
        let c = IntersectionMatrix::COVERED;
        assert!(c.covers() && c.intersects() && c.interiors_meet() && !c.touches() && !c.covered_by());
        let d = IntersectionMatrix::DISJOINT;
        assert!(!d.intersects() && !d.covers() && !d.touches());
        let t = IntersectionMatrix(*b"FF2F01212");
        assert!(t.touches() && t.matches("F***T****") && !t.matches("T********"));
        assert_eq!(c.to_string(), "212FF1FF2");
        assert!(!c.matches("212"));
    }

    #[test]
    fn run_index_merges_overlaps() {
        use crate::output::GridRun;
        let mut r = BurnResult::default();
        r.runs.push(GridRun { row: 2, col_start: 0, col_end: 10, id: 1 });
        r.runs.push(GridRun { row: 2, col_start: 3, col_end: 4, id: 1 });
        let ix = RunIndex::new(&r);
        assert!(ix.covered(1, 2, 5) && ix.covered(1, 2, 0) && !ix.covered(1, 2, 10) && !ix.covered(0, 2, 5));
    }
}
