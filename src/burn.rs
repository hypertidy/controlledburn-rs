//! Public entry points: `burn` and `burn_wkb`.

use std::borrow::Cow;
use std::ops::Range;

use crate::approx::process_polygon_approx;
use crate::coalesce::coalesce_sorted_from;
use crate::ee::bbox::BBox;
use crate::error::BurnError;
use crate::geometry::{Geometry, Polygon};
use crate::grid::{Grid, GridSpec};
use crate::line::process_line;
use crate::output::{BurnResult, GridRun};
use crate::point::process_point;
use crate::polygon::{process_polygon, CoverageScratch};
use crate::wkb::parse_wkb;

/// How polygon boundary cells are classified.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum BurnMode {
    /// Exact coverage fractions via analytical traversal. Boundary cells
    /// appear in `edges` with a fraction in (0, 1); total burned coverage
    /// equals the exact polygon area. The default.
    #[default]
    Coverage,
    /// Cell-centre rule (fasterize semantics): a boundary cell is inside
    /// iff its centre is inside the polygon (left-inclusive). Inside
    /// cells become runs, outside cells are dropped, no `edges`.
    Approx,
}

/// Options for a burn. Build one from [`BurnOptions::default`],
/// [`coverage`](BurnOptions::coverage), [`approx`](BurnOptions::approx) or
/// [`parity`](BurnOptions::parity) and the `with_*` methods. The fields
/// can be read and set, but the struct is non-exhaustive, so it can't be
/// written as a literal outside this crate: a later version may add
/// options without breaking callers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub struct BurnOptions {
    pub mode: BurnMode,
    /// Merge adjacent runs on a row (default `true`). With `false` the
    /// `runs` table is record-for-record identical to the C++ core, which
    /// emits one run per fully covered boundary cell.
    pub coalesce_runs: bool,
    /// One record per cell for a multi-part geometry (default `true`).
    /// A cell that several parts cross gets one edge with the parts'
    /// fractions summed, or a run once the sum reaches 1 (two halves of
    /// a cell either side of a shared border make a full cell); an edge
    /// in a cell another part covers fully is dropped; the geometry's
    /// runs are sorted and coalesced. A cell several line parts cross gets
    /// one record with their lengths summed. Overlapping parts are not
    /// unioned: their sum is capped at a full cell. With `false` each part
    /// keeps its own records, as the C++ core emits them.
    #[cfg_attr(feature = "serde", serde(default = "yes"))]
    pub merge_parts: bool,
    /// Added to each geometry's position to give its `id` (default 0), so
    /// that burns of consecutive slices of a collection, say one per
    /// worker, concatenate with no renumbering. Must be non-negative, and
    /// the last id must fit in an `i32`.
    #[cfg_attr(feature = "serde", serde(default))]
    pub id_offset: i32,
}

#[cfg(feature = "serde")]
fn yes() -> bool {
    true
}

impl Default for BurnOptions {
    fn default() -> Self {
        BurnOptions { mode: BurnMode::Coverage, coalesce_runs: true, merge_parts: true, id_offset: 0 }
    }
}

impl BurnOptions {
    /// Coverage mode with coalescing on.
    pub fn coverage() -> Self {
        Self::default()
    }

    /// Approx mode with coalescing on.
    pub fn approx() -> Self {
        BurnOptions { mode: BurnMode::Approx, ..Self::default() }
    }

    /// Record-for-record parity with the C++ core: coalescing and the
    /// merging of parts off.
    pub fn parity(mode: BurnMode) -> Self {
        BurnOptions { mode, coalesce_runs: false, merge_parts: false, id_offset: 0 }
    }

    pub fn with_mode(mut self, mode: BurnMode) -> Self {
        self.mode = mode;
        self
    }

    pub fn with_coalesce_runs(mut self, on: bool) -> Self {
        self.coalesce_runs = on;
        self
    }

    pub fn with_merge_parts(mut self, on: bool) -> Self {
        self.merge_parts = on;
        self
    }

    pub fn with_id_offset(mut self, offset: i32) -> Self {
        self.id_offset = offset;
        self
    }

    /// The id of the geometry at `k`.
    fn id(&self, k: usize) -> Result<i32, BurnError> {
        i32::try_from(k).ok().and_then(|k| self.id_offset.checked_add(k)).ok_or_else(|| {
            BurnError::InvalidOptions(format!("geometry {k} with id_offset {}: id beyond i32::MAX", self.id_offset))
        })
    }

    fn validate(&self) -> Result<(), BurnError> {
        if self.id_offset < 0 {
            return Err(BurnError::InvalidOptions(format!("id_offset {} is negative", self.id_offset)));
        }
        Ok(())
    }
}

impl From<BurnMode> for BurnOptions {
    fn from(mode: BurnMode) -> Self {
        BurnOptions::default().with_mode(mode)
    }
}

fn ring_needs_closing(ring: &[crate::geometry::Coord]) -> bool {
    ring.len() >= 2 && ring[0] != ring[ring.len() - 1]
}

fn close_polygon(p: &Polygon) -> Polygon {
    Polygon {
        rings: p
            .rings
            .iter()
            .map(|r| {
                let mut r = r.clone();
                if ring_needs_closing(&r) {
                    r.push(r[0]);
                }
                r
            })
            .collect(),
    }
}

/// Close any unclosed polygon ring by appending its first coordinate.
/// Returns a borrowed geometry when nothing needs changing.
fn normalize(g: &Geometry) -> Cow<'_, Geometry> {
    let needs = g.polygons().iter().any(|p| p.rings.iter().any(|r| ring_needs_closing(r)));
    if !needs {
        return Cow::Borrowed(g);
    }
    Cow::Owned(match g {
        Geometry::Polygon(p) => Geometry::Polygon(close_polygon(p)),
        Geometry::MultiPolygon(ps) => Geometry::MultiPolygon(ps.iter().map(close_polygon).collect()),
        other => other.clone(),
    })
}

/// A window of the grid: rows `r0..r1`, columns `c0..c1`, and `near`,
/// its box widened by two cells each side (lines are walked only there).
struct Win {
    r0: i32,
    r1: i32,
    c0: i32,
    c1: i32,
    rows: Range<usize>,
    cols: Range<usize>,
    near: BBox,
}

impl Win {
    fn new(gs: &GridSpec, rows: &Range<u32>, cols: &Range<u32>) -> Self {
        let (dx, dy) = (gs.dx(), gs.dy());
        let x = |c: f64| gs.xmin + c * dx;
        let y = |r: f64| gs.ymax - r * dy;
        let (r0, r1) = (f64::from(rows.start) - 2.0, f64::from(rows.end) + 2.0);
        let (c0, c1) = (f64::from(cols.start) - 2.0, f64::from(cols.end) + 2.0);
        Win {
            near: BBox::new(x(c0), y(r1), x(c1), y(r0)),
            r0: rows.start as i32,
            r1: rows.end as i32,
            c0: cols.start as i32,
            c1: cols.end as i32,
            rows: rows.start as usize..rows.end as usize,
            cols: cols.start as usize..cols.end as usize,
        }
    }

    fn has(&self, row: i32, col: i32) -> bool {
        (self.r0..self.r1).contains(&row) && (self.c0..self.c1).contains(&col)
    }
}

/// Keep the records from `run_start`, `edge_start` and `line_start` on
/// that lie in `win`, cutting runs at its columns.
fn clip(out: &mut BurnResult, win: &Win, run_start: usize, edge_start: usize, line_start: usize) {
    let mut w = run_start;
    for i in run_start..out.runs.len() {
        let mut r = out.runs[i];
        r.col_start = r.col_start.max(win.c0);
        r.col_end = r.col_end.min(win.c1);
        if (win.r0..win.r1).contains(&r.row) && r.col_start < r.col_end {
            out.runs[w] = r;
            w += 1;
        }
    }
    out.runs.truncate(w);
    let mut i = edge_start;
    out.edges.retain(|e| {
        i += 1;
        i <= edge_start || win.has(e.row, e.col)
    });
    let mut i = line_start;
    out.lines.retain(|l| {
        i += 1;
        i <= line_start || win.has(l.row, l.col)
    });
}

#[allow(clippy::too_many_arguments)]
fn process_geometry(
    g: &Geometry,
    full_grid: &Grid<0>,
    gs: &GridSpec,
    geom_id: i32,
    out: &mut BurnResult,
    opts: &BurnOptions,
    sc: &mut CoverageScratch,
    win: Option<&Win>,
) -> Result<(), String> {
    let run_start = out.runs.len();
    let edge_start = out.edges.len();
    let line_start = out.lines.len();
    let polys = g.polygons();
    note_short_rings(polys, geom_id, out);
    for poly in polys {
        match opts.mode {
            BurnMode::Approx => {
                let rows = win.map(|w| (i64::from(w.r0), i64::from(w.r1) - 1));
                process_polygon_approx(poly, gs, geom_id, &mut out.runs, rows)
            }
            BurnMode::Coverage => {
                let w = win.map(|w| (&w.rows, &w.cols));
                process_polygon(poly, full_grid, geom_id, &mut out.runs, &mut out.edges, sc, w)?
            }
        }
    }
    if let Some(w) = win {
        clip(out, w, run_start, edge_start, line_start);
    }
    let merge = opts.merge_parts && polys.len() > 1;
    if merge {
        merge_part_edges(out, edge_start);
    }
    if (opts.coalesce_runs || merge) && out.runs.len() > run_start {
        // One polygon's runs come out of the sweep ordered by (row,
        // col_start); a multipolygon's are in that order only per part, so
        // sort them together first, or a later part's run left of an
        // earlier part's on a shared row is swallowed by the merge.
        if polys.len() > 1 {
            out.runs[run_start..].sort_by_key(|r| (r.row, r.col_start));
        }
        coalesce_sorted_from(&mut out.runs, run_start);
    }
    if merge {
        drop_edges_in_runs(out, run_start, edge_start);
    }
    for (i, line) in g.lines().iter().enumerate() {
        if line.len() == 1 {
            out.note(geom_id, format!("skipped line part {i}: 1 coordinate, 2 needed"));
        }
        process_line(line, full_grid, geom_id, &mut out.lines, win.map(|w| &w.near))?;
    }
    if let Some(w) = win {
        clip(out, w, out.runs.len(), out.edges.len(), line_start);
    }
    if opts.merge_parts && g.lines().len() > 1 {
        merge_part_lines(out, line_start);
    }
    let point_start = out.points.len();
    for pt in g.points() {
        process_point(pt, full_grid, geom_id, &mut out.points);
    }
    if let Some(w) = win {
        let mut i = point_start;
        out.points.retain(|p| {
            i += 1;
            i <= point_start || w.has(p.row, p.col)
        });
    }
    Ok(())
}

/// A note for each ring the engines skip as too short: fewer than 4
/// coordinates once closed. A short exterior skips its part, holes and
/// all; a short hole skips only itself. Empty rings are skipped silently.
fn note_short_rings(polys: &[Polygon], id: i32, out: &mut BurnResult) {
    for (pi, p) in polys.iter().enumerate() {
        for (ri, r) in p.rings.iter().enumerate() {
            if r.is_empty() || r.len() >= 4 {
                continue;
            }
            let n = r.len();
            if ri == 0 {
                out.note(
                    id,
                    format!("skipped polygon part {pi}: exterior ring has {n} coordinates when closed, 4 needed"),
                );
                break;
            }
            out.note(id, format!("skipped hole {ri} of polygon part {pi}: {n} coordinates when closed, 4 needed"));
        }
    }
}

/// The sweep's tolerance on a cell's coverage (`process_polygon`).
const TOL: f32 = 1e-6;

/// One edge per cell for the edges of one geometry, `edges[start..]`:
/// sum the parts' fractions per cell, and make a cell whose sum reaches
/// 1 a one-cell run.
fn merge_part_edges(out: &mut BurnResult, start: usize) {
    let edges = &mut out.edges;
    // stable: a cell's fractions are summed in part order
    edges[start..].sort_by_key(|e| (e.row, e.col));
    let mut w = start;
    for i in start..edges.len() {
        let cur = edges[i];
        if w > start {
            let last = &mut edges[w - 1];
            if last.row == cur.row && last.col == cur.col {
                last.fraction += cur.fraction;
                continue;
            }
        }
        edges[w] = cur;
        w += 1;
    }
    edges.truncate(w);
    let mut w = start;
    for i in start..edges.len() {
        let e = edges[i];
        if e.fraction >= 1.0 - TOL {
            out.runs.push(GridRun { row: e.row, col_start: e.col, col_end: e.col + 1, id: e.id });
        } else {
            edges[w] = e;
            w += 1;
        }
    }
    edges.truncate(w);
}

/// Drop the edges of one geometry (`edges[edge_start..]`, sorted by
/// cell) that lie in one of its runs (`runs[run_start..]`, sorted and
/// coalesced): another part covers the cell fully.
fn drop_edges_in_runs(out: &mut BurnResult, run_start: usize, edge_start: usize) {
    let runs = &out.runs[run_start..];
    let edges = &mut out.edges;
    let mut j = 0;
    let mut w = edge_start;
    for i in edge_start..edges.len() {
        let e = edges[i];
        while j < runs.len() && (runs[j].row, runs[j].col_end) <= (e.row, e.col) {
            j += 1;
        }
        let covered = j < runs.len() && runs[j].row == e.row && runs[j].col_start <= e.col;
        if !covered {
            edges[w] = e;
            w += 1;
        }
    }
    edges.truncate(w);
}

/// One record per cell for the line parts of one geometry,
/// `lines[start..]`: lengths summed per cell, in part order.
fn merge_part_lines(out: &mut BurnResult, start: usize) {
    let lines = &mut out.lines;
    lines[start..].sort_by_key(|l| (l.row, l.col));
    let mut w = start;
    for i in start..lines.len() {
        let cur = lines[i];
        if w > start {
            let last = &mut lines[w - 1];
            if last.row == cur.row && last.col == cur.col {
                last.length += cur.length;
                continue;
            }
        }
        lines[w] = cur;
        w += 1;
    }
    lines.truncate(w);
}

#[allow(clippy::too_many_arguments)]
fn burn_one(
    g: &Geometry,
    id: i32,
    full_grid: &Grid<0>,
    gs: &GridSpec,
    out: &mut BurnResult,
    opts: &BurnOptions,
    sc: &mut CoverageScratch,
    win: Option<&Win>,
) {
    if g.is_empty() {
        return;
    }
    if !g.is_finite() {
        out.note(id, "skipped geometry with non-finite coordinates");
        return;
    }
    let g = normalize(g);
    if let Err(e) = process_geometry(&g, full_grid, gs, id, out, opts, sc, win) {
        out.note(id, format!("error processing geometry: {e}"));
    }
}

/// Burn a set of geometries onto the grid. Geometry `k` (its 0-based
/// position in `geoms`) is `id = k + opts.id_offset` in the output tables
/// and in `notes`.
///
/// Empty geometries are skipped silently (an all-NaN point, WKB's
/// `POINT EMPTY`, is empty). Geometries with non-finite coordinates and
/// geometries the engine cannot process are skipped with a note; so are
/// rings of fewer than 4 coordinates (once closed) and lines of one
/// coordinate, the rest of the geometry being burned. Unclosed polygon
/// rings are closed. The only errors are an invalid `GridSpec` and an
/// invalid `id_offset`.
pub fn burn(geoms: &[Geometry], grid: &GridSpec, opts: impl Into<BurnOptions>) -> Result<BurnResult, BurnError> {
    grid.validate()?;
    let opts = opts.into();
    opts.validate()?;
    if let Some(last) = geoms.len().checked_sub(1) {
        opts.id(last)?;
    }
    let full_grid = grid.bounded();
    let mut out = BurnResult::default();
    let mut sc = CoverageScratch::default();
    for (k, g) in geoms.iter().enumerate() {
        burn_one(g, opts.id(k)?, &full_grid, grid, &mut out, &opts, &mut sc, None);
    }
    Ok(out)
}

/// Burn geometries supplied as WKB blobs (one blob per geometry; ISO and
/// EWKB, both byte orders, Z/M skipped). Empty blobs are skipped
/// silently, unparseable blobs with a note.
pub fn burn_wkb<'a, I>(wkb: I, grid: &GridSpec, opts: impl Into<BurnOptions>) -> Result<BurnResult, BurnError>
where
    I: IntoIterator<Item = &'a [u8]>,
{
    grid.validate()?;
    let opts = opts.into();
    opts.validate()?;
    let full_grid = grid.bounded();
    let mut out = BurnResult::default();
    let mut sc = CoverageScratch::default();
    for (k, blob) in wkb.into_iter().enumerate() {
        if blob.is_empty() {
            continue;
        }
        let id = opts.id(k)?;
        match parse_wkb(blob) {
            Ok(g) => burn_one(&g, id, &full_grid, grid, &mut out, &opts, &mut sc, None),
            Err(e) => out.note(id, format!("failed to parse WKB: {e}")),
        }
    }
    Ok(out)
}

/// Burn `geoms` onto a window of `grid`: rows `rows` and columns `cols`.
/// The records are those [`burn`] gives in the window's cells, with the
/// same values, in the same order, and with the full grid's indices: runs
/// are cut at the window's columns, and nothing outside it is kept. The
/// windows of a tiling put together are the whole-grid burn, so each tile
/// can burn just the geometries that meet it, on its own thread or
/// machine (see `design/tile-parallel.md`).
///
/// The window keeps the full grid's lattice and edge rules: every wall
/// and cell centre is the full grid's, a line on a wall stays in the cell
/// the whole-grid burn gives it, and a point on a seam between windows is
/// in exactly one of them (only the grid's own right and bottom edges are
/// closed). [`slice_tiles`](crate::slice_tiles) of the whole burn gives
/// the same records, indexed on each tile.
///
/// Notes are as [`burn`] gives them, for every geometry in `geoms`.
///
/// Each geometry is burned as [`burn`] burns it, on the same sub-grid, so
/// the values are the same computation, not a recomputation: a polygon's
/// rows outside the window are skipped (they are independent) but its
/// boundary is walked whole; a line is walked only along its runs of
/// segments that pass within two cells of the window (the walker computes
/// each crossing from the segment's own vertices, so a walk started
/// mid-line agrees with the whole walk away from its start); points and
/// the columns outside the window are filtered. So a window costs about
/// the vertices of the polygons that meet it, plus its own output. For a
/// polygon much larger than the window, burning once and slicing is
/// cheaper.
///
/// The errors are those of [`burn`], and an empty window or one outside
/// the grid (`InvalidGrid`).
///
/// # Examples
///
/// ```
/// use controlledburn::{burn, burn_window, BurnOptions, Coord, Geometry,
///     GridSpec, Polygon};
///
/// let square = Geometry::Polygon(Polygon::new(vec![vec![
///     Coord::new(2.5, 4.5), Coord::new(6.5, 4.5), Coord::new(6.5, 8.5),
///     Coord::new(2.5, 8.5), Coord::new(2.5, 4.5),
/// ]]));
/// let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
/// let geoms = [square];
///
/// // the left and right halves of the grid
/// let left = burn_window(&geoms, &grid, 0..10, 0..5, BurnOptions::default()).unwrap();
/// let right = burn_window(&geoms, &grid, 0..10, 5..10, BurnOptions::default()).unwrap();
/// let whole = burn(&geoms, &grid, BurnOptions::default()).unwrap();
/// assert!(left.runs.iter().all(|r| r.col_end <= 5));
/// assert_eq!(left.covered_cells() + right.covered_cells(), whole.covered_cells());
/// ```
pub fn burn_window(
    geoms: &[Geometry],
    grid: &GridSpec,
    rows: Range<u32>,
    cols: Range<u32>,
    opts: impl Into<BurnOptions>,
) -> Result<BurnResult, BurnError> {
    grid.validate()?;
    if rows.is_empty() || cols.is_empty() || rows.end > grid.nrow || cols.end > grid.ncol {
        return Err(BurnError::InvalidGrid(format!(
            "window rows {rows:?}, cols {cols:?} must be non-empty and within {} rows, {} cols",
            grid.nrow, grid.ncol
        )));
    }
    let opts = opts.into();
    opts.validate()?;
    if let Some(last) = geoms.len().checked_sub(1) {
        opts.id(last)?;
    }
    let win = Win::new(grid, &rows, &cols);
    let full_grid = grid.bounded();
    let mut out = BurnResult::default();
    let mut sc = CoverageScratch::default();
    for (k, g) in geoms.iter().enumerate() {
        burn_one(g, opts.id(k)?, &full_grid, grid, &mut out, &opts, &mut sc, Some(&win));
    }
    Ok(out)
}

/// Burn a sequence of [`geo-traits`](https://docs.rs/geo-traits)
/// geometries (feature `geo-traits`): geo-types, the `wkb` crate's reader,
/// GeoArrow arrays, or anything else that implements `GeometryTrait`.
/// Geometry `k` (its 0-based position) is `id = k + opts.id_offset`, as
/// in [`burn`].
///
/// Each geometry is copied with
/// [`Geometry::from_geo_traits`](crate::Geometry::from_geo_traits) and
/// burned as [`burn`] burns it, so the output is the same as converting
/// first. A `GeometryCollection` is skipped with a note.
///
/// ```
/// use controlledburn::{burn_geo_traits, BurnOptions, GridSpec};
///
/// let ring = vec![(2.5, 4.5), (6.5, 4.5), (6.5, 8.5), (2.5, 8.5)];
/// let square = geo_types::Polygon::new(ring.into(), vec![]);
/// let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
/// let r = burn_geo_traits([&square], &grid, BurnOptions::default()).unwrap();
/// assert!((r.covered_cells() - 16.0).abs() < 1e-5);
/// ```
#[cfg(feature = "geo-traits")]
pub fn burn_geo_traits<I, G>(geoms: I, grid: &GridSpec, opts: impl Into<BurnOptions>) -> Result<BurnResult, BurnError>
where
    I: IntoIterator<Item = G>,
    G: geo_traits::GeometryTrait<T = f64>,
{
    grid.validate()?;
    let opts = opts.into();
    opts.validate()?;
    let full_grid = grid.bounded();
    let mut out = BurnResult::default();
    let mut sc = CoverageScratch::default();
    for (k, g) in geoms.into_iter().enumerate() {
        let id = opts.id(k)?;
        match Geometry::from_geo_traits(&g) {
            Ok(g) => burn_one(&g, id, &full_grid, grid, &mut out, &opts, &mut sc, None),
            Err(e) => out.note(id, format!("skipped geometry: {e}")),
        }
    }
    Ok(out)
}
