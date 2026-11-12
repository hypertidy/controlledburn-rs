//! Public entry points: `burn` and `burn_wkb`.

use std::borrow::Cow;
use std::ops::Range;

use crate::approx::process_polygon_approx;
use crate::coalesce::coalesce_sorted_from;
use crate::ee::bbox::BBox;
use crate::error::BurnError;
use crate::geometry::{Coord, Geometry, Polygon};
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
    pub(crate) fn id(&self, k: usize) -> Result<i32, BurnError> {
        i32::try_from(k).ok().and_then(|k| self.id_offset.checked_add(k)).ok_or_else(|| {
            BurnError::InvalidOptions(format!("geometry {k} with id_offset {}: id beyond i32::MAX", self.id_offset))
        })
    }

    pub(crate) fn validate(&self) -> Result<(), BurnError> {
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

/// `burn_par` uses at most one thread per this much work (a geometry's
/// work being its coordinates plus the grid rows it spans), and burns on
/// the calling thread below two.
const WORK_PER_THREAD: usize = 1 << 13;

/// Chunks of geometries per thread: enough that a costly chunk taken last
/// does not leave the other threads idle for long.
const CHUNKS_PER_THREAD: usize = 16;

/// [`burn`] on `threads` threads (0: all available). The result is the
/// same as `burn`'s, record for record and in the same order, notes
/// included.
///
/// The geometries are split into contiguous chunks of about equal work
/// (coordinates plus the grid rows each geometry spans), several per
/// thread. Each thread takes the costliest
/// chunk left, burns it into its own tables with its own scratch buffers,
/// and takes the next; the chunks' tables are then put together in
/// geometry order. The extra memory is the chunks' tables, which are freed
/// as they are copied.
///
/// One geometry is never split, so a collection dominated by one large
/// geometry gets little speedup: burn that once, or burn tiles of it with
/// [`burn_window`].
///
/// It uses at most one thread per 8,192 (2^13) of that work, so `threads`
/// is an upper bound; with one thread or one geometry it runs [`burn`] on
/// the calling thread. Threads come from `std::thread::scope`; on a target
/// without threads (such as `wasm32-unknown-unknown`) pass `threads = 1`,
/// or use [`burn`].
///
/// # Examples
///
/// ```
/// use controlledburn::{burn, burn_par, BurnOptions, Coord, Geometry,
///     GridSpec, Polygon};
///
/// let squares: Vec<Geometry> = (0..20)
///     .map(|i| {
///         let (x, y) = ((i % 5) as f64 * 2.0 + 0.3, (i / 5) as f64 * 2.0 + 0.3);
///         Geometry::Polygon(Polygon::new(vec![vec![
///             Coord::new(x, y), Coord::new(x + 1.4, y), Coord::new(x + 1.4, y + 1.4),
///             Coord::new(x, y + 1.4), Coord::new(x, y),
///         ]]))
///     })
///     .collect();
/// let grid = GridSpec::new(0.0, 0.0, 10.0, 8.0, 20, 16);
/// let par = burn_par(&squares, &grid, BurnOptions::default(), 4).unwrap();
/// assert_eq!(par, burn(&squares, &grid, BurnOptions::default()).unwrap());
/// ```
pub fn burn_par(
    geoms: &[Geometry],
    grid: &GridSpec,
    opts: impl Into<BurnOptions>,
    threads: usize,
) -> Result<BurnResult, BurnError> {
    burn_par_from(geoms, grid, opts.into(), threads, WORK_PER_THREAD)
}

/// `burn_par` with one thread per `per_thread` of work at most; 0
/// means no cap (the tests force the parallel path on small inputs).
fn burn_par_from(
    geoms: &[Geometry],
    grid: &GridSpec,
    opts: BurnOptions,
    threads: usize,
    per_thread: usize,
) -> Result<BurnResult, BurnError> {
    let cost: Vec<usize> = geoms.iter().map(|g| work(g, grid)).collect();
    let Some(threads) = par_threads(&cost, threads, per_thread) else {
        return burn(geoms, grid, opts);
    };
    grid.validate()?;
    opts.validate()?;
    opts.id(geoms.len() - 1)?;
    Ok(par_chunks(&cost, grid, &opts, threads, |k, id, full_grid, out, sc| {
        burn_one(&geoms[k], id, full_grid, grid, out, &opts, sc, None)
    }))
}

/// [`burn_wkb`] on `threads` threads (0: all available): [`burn_par`] for
/// WKB blobs, with the same result as `burn_wkb`, record for record and in
/// the same order, notes included.
///
/// Each blob is parsed by the thread that burns it, as it burns it, so no
/// parsed geometry outlives its burn and none crosses a thread. Parsing
/// first and calling [`burn_par`] instead allocates every geometry at once
/// and frees it afterwards, which for many small geometries can cost as
/// much as the burn. The work of each blob, for balancing the threads, is
/// read from it without building the geometry.
///
/// Chunking, the thread cap and the fallback to the calling thread are
/// those of [`burn_par`]. The errors are those of [`burn_wkb`], with ids
/// checked for the whole collection first.
///
/// # Examples
///
/// ```
/// use controlledburn::{burn_wkb, burn_wkb_par, BurnOptions, GridSpec};
///
/// // little-endian WKB polygons: small squares
/// let square = |x: f64, y: f64| {
///     let mut b = vec![1u8];
///     b.extend(3u32.to_le_bytes());
///     b.extend(1u32.to_le_bytes());
///     b.extend(5u32.to_le_bytes());
///     for (px, py) in [(x, y), (x + 1.4, y), (x + 1.4, y + 1.4), (x, y + 1.4), (x, y)] {
///         b.extend(px.to_le_bytes());
///         b.extend(py.to_le_bytes());
///     }
///     b
/// };
/// let blobs: Vec<Vec<u8>> = (0..20).map(|i| square((i % 5) as f64 * 2.0 + 0.3, (i / 5) as f64 * 2.0 + 0.3)).collect();
/// let grid = GridSpec::new(0.0, 0.0, 10.0, 8.0, 20, 16);
/// let par = burn_wkb_par(&blobs, &grid, BurnOptions::default(), 4).unwrap();
/// let serial = burn_wkb(blobs.iter().map(|b| b.as_slice()), &grid, BurnOptions::default()).unwrap();
/// assert_eq!(par, serial);
/// ```
pub fn burn_wkb_par<B: AsRef<[u8]> + Sync>(
    wkb: &[B],
    grid: &GridSpec,
    opts: impl Into<BurnOptions>,
    threads: usize,
) -> Result<BurnResult, BurnError> {
    burn_wkb_par_from(wkb, grid, opts.into(), threads, WORK_PER_THREAD)
}

fn burn_wkb_par_from<B: AsRef<[u8]> + Sync>(
    wkb: &[B],
    grid: &GridSpec,
    opts: BurnOptions,
    threads: usize,
    per_thread: usize,
) -> Result<BurnResult, BurnError> {
    let cost: Vec<usize> = wkb.iter().map(|b| work_wkb(b.as_ref(), grid)).collect();
    let Some(threads) = par_threads(&cost, threads, per_thread) else {
        return burn_wkb(wkb.iter().map(|b| b.as_ref()), grid, opts);
    };
    grid.validate()?;
    opts.validate()?;
    opts.id(wkb.len() - 1)?;
    Ok(par_chunks(&cost, grid, &opts, threads, |k, id, full_grid, out, sc| {
        let blob = wkb[k].as_ref();
        if blob.is_empty() {
            return;
        }
        match parse_wkb(blob) {
            Ok(g) => burn_one(&g, id, full_grid, grid, out, &opts, sc, None),
            Err(e) => out.note(id, format!("failed to parse WKB: {e}")),
        }
    }))
}

/// The threads to use for items of `cost`, or `None` to burn on the
/// calling thread (fewer than two threads, or fewer than two items).
fn par_threads(cost: &[usize], threads: usize, per_thread: usize) -> Option<usize> {
    let threads = if threads == 0 { std::thread::available_parallelism().map_or(1, |n| n.get()) } else { threads };
    let total: usize = cost.iter().sum();
    let threads = total.checked_div(per_thread).map_or(threads, |n| threads.min(n));
    (threads >= 2 && cost.len() >= 2).then_some(threads)
}

/// Burn items `0..cost.len()` on `threads` threads in contiguous chunks of
/// about equal cost, costliest first, and join the chunks' tables in item
/// order. `burn_k(k, id, full_grid, out, scratch)` burns item `k`; ids
/// must have been checked.
fn par_chunks<F>(cost: &[usize], grid: &GridSpec, opts: &BurnOptions, threads: usize, burn_k: F) -> BurnResult
where
    F: Fn(usize, i32, &Grid<0>, &mut BurnResult, &mut CoverageScratch) + Sync,
{
    let total: usize = cost.iter().sum();
    let target = (total + threads * CHUNKS_PER_THREAD - 1) / (threads * CHUNKS_PER_THREAD);
    let mut chunks: Vec<Range<usize>> = Vec::new();
    let (mut start, mut acc) = (0, 0);
    for (k, &c) in cost.iter().enumerate() {
        acc += c;
        if acc >= target {
            chunks.push(start..k + 1);
            (start, acc) = (k + 1, 0);
        }
    }
    if start < cost.len() {
        chunks.push(start..cost.len());
    }
    // Costliest first, so the last chunks taken are small.
    let chunk_cost = |r: &Range<usize>| cost[r.clone()].iter().sum::<usize>();
    let mut order: Vec<usize> = (0..chunks.len()).collect();
    order.sort_by_key(|&i| std::cmp::Reverse(chunk_cost(&chunks[i])));

    let full_grid = grid.bounded();
    let next = std::sync::atomic::AtomicUsize::new(0);
    let (chunks, order, next, full_grid, burn_k) = (&chunks, &order, &next, &full_grid, &burn_k);
    let mut parts: Vec<(usize, BurnResult)> = std::thread::scope(|sc| {
        let hs: Vec<_> = (0..threads.min(chunks.len()))
            .map(|_| {
                sc.spawn(move || {
                    let mut mine = Vec::new();
                    let mut sc = CoverageScratch::default();
                    loop {
                        let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        let Some(&c) = order.get(i) else { break };
                        let mut out = BurnResult::default();
                        for k in chunks[c].clone() {
                            // ids were checked by the caller
                            let id = opts.id(k).expect("id checked");
                            burn_k(k, id, full_grid, &mut out, &mut sc);
                        }
                        mine.push((c, out));
                    }
                    mine
                })
            })
            .collect();
        hs.into_iter().flat_map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e))).collect()
    });
    parts.sort_unstable_by_key(|p| p.0);

    let mut out = BurnResult::default();
    out.runs.reserve_exact(parts.iter().map(|p| p.1.runs.len()).sum());
    out.edges.reserve_exact(parts.iter().map(|p| p.1.edges.len()).sum());
    out.lines.reserve_exact(parts.iter().map(|p| p.1.lines.len()).sum());
    out.points.reserve_exact(parts.iter().map(|p| p.1.points.len()).sum());
    for (_, p) in parts {
        out.runs.extend_from_slice(&p.runs);
        out.edges.extend_from_slice(&p.edges);
        out.lines.extend_from_slice(&p.lines);
        out.points.extend_from_slice(&p.points);
        out.notes.extend(p.notes);
    }
    out
}

/// A geometry's estimated burn cost: one, plus its coordinates (the
/// walk), plus the grid rows its extent spans (the sweep and the output).
/// The extent is taken from at most 16 coordinates of each part, so the
/// estimate costs O(parts), not O(coordinates).
fn work(g: &Geometry, grid: &GridSpec) -> usize {
    let (mut n, mut lo, mut hi) = (0usize, f64::INFINITY, f64::NEG_INFINITY);
    let mut part = |c: &[Coord]| {
        n += c.len();
        let step = (c.len() + 15) / 16 + usize::from(c.is_empty());
        for p in c.iter().step_by(step).chain(c.last()) {
            lo = lo.min(p.y);
            hi = hi.max(p.y);
        }
    };
    part(g.points());
    g.lines().iter().for_each(|l| part(l));
    g.polygons().iter().filter_map(|p| p.rings.first()).for_each(|r| part(r));
    n += g.polygons().iter().flat_map(|p| p.rings.iter().skip(1)).map(Vec::len).sum::<usize>();
    work_from(n, lo, hi, grid)
}

/// [`work`] of a WKB blob, read without building the geometry; 1 for a
/// blob that does not parse.
fn work_wkb(blob: &[u8], grid: &GridSpec) -> usize {
    crate::wkb::sample_extent(blob).map_or(1, |(n, lo, hi)| work_from(n, lo, hi, grid))
}

fn work_from(n: usize, lo: f64, hi: f64, grid: &GridSpec) -> usize {
    let (lo, hi) = (lo.max(grid.ymin), hi.min(grid.ymax));
    let rows = if hi > lo { ((hi - lo) / grid.dy()).min(grid.nrow as f64) as usize + 1 } else { 0 };
    1 + n + rows
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Polygon;
    use proptest::prelude::*;

    fn ring(cx: f64, cy: f64, r: f64, k: usize, wobble: &[f64]) -> Vec<Coord> {
        let mut v: Vec<Coord> = (0..k)
            .map(|i| {
                let a = i as f64 / k as f64 * std::f64::consts::TAU;
                let rr = r * (0.6 + 0.4 * wobble[i % wobble.len()]);
                Coord::new(cx + rr * a.cos(), cy + rr * a.sin())
            })
            .collect();
        v.push(v[0]);
        v
    }

    proptest! {
        /// The parallel burn equals the serial one, record for record and
        /// note for note, whatever the threads, options and input (the
        /// parallel path forced on from two geometries).
        #[test]
        fn parallel_equals_serial(
            centres in prop::collection::vec((-3.0..40.0f64, -3.0..26.0f64, 0.3..6.0f64, 0u8..6), 0..40),
            wobble in prop::collection::vec(0.0..1.0f64, 9),
            threads in 2usize..12, approx in any::<bool>(), coalesce in any::<bool>(),
            merge in any::<bool>(), id_offset in 0i32..1000,
        ) {
            let grid = GridSpec::new(0.0, 0.0, 37.0, 23.0, 53, 31);
            let geoms: Vec<Geometry> = centres
                .iter()
                .map(|&(x, y, r, kind)| match kind {
                    0 => Geometry::Polygon(Polygon::new(vec![ring(x, y, r, 9, &wobble)])),
                    1 => Geometry::MultiPolygon(vec![
                        Polygon::new(vec![ring(x, y, r, 7, &wobble), ring(x, y, r * 0.3, 5, &wobble)]),
                        Polygon::new(vec![ring(x + r * 0.8, y, r * 0.7, 6, &wobble)]),
                    ]),
                    2 => Geometry::LineString(ring(x, y, r, 8, &wobble)),
                    3 => Geometry::MultiPoint((0..7).map(|i| Coord::new(x + i as f64 * 0.7, y - i as f64 * 0.4)).collect()),
                    4 => Geometry::Polygon(Polygon::new(vec![])),
                    _ => Geometry::Point(Coord::new(f64::NAN, y)),
                })
                .collect();
            let opts = BurnOptions::default()
                .with_mode(if approx { BurnMode::Approx } else { BurnMode::Coverage })
                .with_coalesce_runs(coalesce)
                .with_merge_parts(merge)
                .with_id_offset(id_offset);
            let want = burn(&geoms, &grid, opts).unwrap();
            let got = burn_par_from(&geoms, &grid, opts, threads, 0).unwrap();
            prop_assert_eq!(got, want);
        }
    }

    /// ISO WKB of `g`, in either byte order, with `extra` (0-2) zero
    /// ordinates after x and y.
    fn to_wkb(g: &Geometry, big: bool, extra: usize) -> Vec<u8> {
        fn head(b: &mut Vec<u8>, big: bool, base: u32, extra: usize) {
            b.push(u8::from(!big));
            let code = base + [0, 1000, 3000][extra];
            b.extend(if big { code.to_be_bytes() } else { code.to_le_bytes() });
        }
        fn n(b: &mut Vec<u8>, big: bool, v: usize) {
            b.extend(if big { (v as u32).to_be_bytes() } else { (v as u32).to_le_bytes() });
        }
        fn pt(b: &mut Vec<u8>, big: bool, c: &Coord, extra: usize) {
            for v in [c.x, c.y].into_iter().chain(std::iter::repeat(0.0).take(extra)) {
                b.extend(if big { v.to_be_bytes() } else { v.to_le_bytes() });
            }
        }
        fn seq(b: &mut Vec<u8>, big: bool, s: &[Coord], extra: usize) {
            n(b, big, s.len());
            s.iter().for_each(|c| pt(b, big, c, extra));
        }
        fn poly(b: &mut Vec<u8>, big: bool, p: &Polygon, extra: usize) {
            n(b, big, p.rings.len());
            p.rings.iter().for_each(|r| seq(b, big, r, extra));
        }
        let mut b = Vec::new();
        match g {
            Geometry::Point(c) => {
                head(&mut b, big, 1, extra);
                pt(&mut b, big, c, extra);
            }
            Geometry::LineString(l) => {
                head(&mut b, big, 2, extra);
                seq(&mut b, big, l, extra);
            }
            Geometry::Polygon(p) => {
                head(&mut b, big, 3, extra);
                poly(&mut b, big, p, extra);
            }
            Geometry::MultiPoint(ps) => {
                head(&mut b, big, 4, extra);
                n(&mut b, big, ps.len());
                for c in ps {
                    head(&mut b, big, 1, extra);
                    pt(&mut b, big, c, extra);
                }
            }
            Geometry::MultiLineString(ls) => {
                head(&mut b, big, 5, extra);
                n(&mut b, big, ls.len());
                for l in ls {
                    head(&mut b, big, 2, extra);
                    seq(&mut b, big, l, extra);
                }
            }
            Geometry::MultiPolygon(ps) => {
                head(&mut b, big, 6, extra);
                n(&mut b, big, ps.len());
                for p in ps {
                    head(&mut b, big, 3, extra);
                    poly(&mut b, big, p, extra);
                }
            }
        }
        b
    }

    /// Geometries of every kind for the WKB tests, with a ring of 40 or
    /// more coordinates (so sampling skips some).
    fn mixed(centres: &[(f64, f64, f64, u8)], wobble: &[f64]) -> Vec<Geometry> {
        centres
            .iter()
            .map(|&(x, y, r, kind)| match kind {
                0 => Geometry::Polygon(Polygon::new(vec![ring(x, y, r, 9, wobble)])),
                1 => Geometry::MultiPolygon(vec![
                    Polygon::new(vec![ring(x, y, r, 41, wobble), ring(x, y, r * 0.3, 5, wobble)]),
                    Polygon::new(vec![ring(x + r * 0.8, y, r * 0.7, 6, wobble)]),
                ]),
                2 => Geometry::LineString(ring(x, y, r, 37, wobble)),
                3 => {
                    Geometry::MultiPoint((0..23).map(|i| Coord::new(x + i as f64 * 0.7, y - i as f64 * 0.4)).collect())
                }
                // an open line whose last coordinate alone reaches higher:
                // 38 coordinates, so the sample's step (3) misses it
                4 => Geometry::MultiLineString(vec![
                    ring(x, y, r, 8, wobble),
                    vec![],
                    (0..38).map(|i| Coord::new(x + i as f64 * 0.2, if i == 37 { y + 5.0 } else { y })).collect(),
                ]),
                5 => Geometry::Polygon(Polygon::new(vec![])),
                6 => Geometry::Point(Coord::new(x, y)),
                _ => Geometry::Point(Coord::new(f64::NAN, y)),
            })
            .collect()
    }

    proptest! {
        /// `burn_wkb_par` equals `burn_wkb`, record for record and note for
        /// note, with blobs in both byte orders, with Z and M, empty and
        /// unparseable blobs among them.
        #[test]
        fn wkb_parallel_equals_serial(
            centres in prop::collection::vec((-3.0..40.0f64, -3.0..26.0f64, 0.3..6.0f64, 0u8..10), 0..40),
            wobble in prop::collection::vec(0.0..1.0f64, 9),
            threads in 2usize..12, approx in any::<bool>(), merge in any::<bool>(), id_offset in 0i32..1000,
        ) {
            let grid = GridSpec::new(0.0, 0.0, 37.0, 23.0, 53, 31);
            let geoms = mixed(&centres, &wobble);
            let blobs: Vec<Vec<u8>> = geoms
                .iter()
                .zip(&centres)
                .enumerate()
                .map(|(k, (g, c))| match c.3 {
                    8 => Vec::new(),
                    9 => vec![1, 3, 0, 0, 0, 7],
                    _ => to_wkb(g, k % 2 == 1, k % 3),
                })
                .collect();
            let opts = BurnOptions::default()
                .with_mode(if approx { BurnMode::Approx } else { BurnMode::Coverage })
                .with_merge_parts(merge)
                .with_id_offset(id_offset);
            let want = burn_wkb(blobs.iter().map(|b| b.as_slice()), &grid, opts).unwrap();
            let got = burn_wkb_par_from(&blobs, &grid, opts, threads, 0).unwrap();
            prop_assert_eq!(got, want);
        }

        /// The work read from a blob is the work of the parsed geometry.
        #[test]
        fn wkb_work_is_work(
            centres in prop::collection::vec((-3.0..40.0f64, -3.0..26.0f64, 0.3..6.0f64, 0u8..8), 1..20),
            wobble in prop::collection::vec(0.0..1.0f64, 9),
        ) {
            let grid = GridSpec::new(0.0, 0.0, 37.0, 23.0, 53, 31);
            for g in mixed(&centres, &wobble) {
                for (big, extra) in [(false, 0), (true, 0), (false, 1), (true, 2)] {
                    prop_assert_eq!(work_wkb(&to_wkb(&g, big, extra), &grid), work(&g, &grid), "{:?}", g);
                }
            }
        }
    }

    #[test]
    fn errors_as_serial() {
        let sq = |x: f64| Geometry::Polygon(Polygon::new(vec![ring(x, 5.0, 3.0, 12, &[0.5])]));
        let geoms = [sq(3.0), sq(6.0)];
        let ok = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
        let bad = GridSpec::new(0.0, 0.0, 10.0, 10.0, 0, 10);
        let cases = [
            (bad, BurnOptions::default()),
            (ok, BurnOptions::default().with_id_offset(-1)),
            (ok, BurnOptions::default().with_id_offset(i32::MAX)),
        ];
        for (grid, opts) in cases {
            let want = burn(&geoms, &grid, opts).unwrap_err().to_string();
            assert_eq!(burn_par_from(&geoms, &grid, opts, 4, 0).unwrap_err().to_string(), want);
        }
    }

    #[test]
    fn falls_back_to_serial() {
        // one thread, one geometry, or few coordinates: same answer
        let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
        let sq = |x: f64| Geometry::Polygon(Polygon::new(vec![ring(x, 5.0, 3.0, 12, &[0.5])]));
        let want = burn(&[sq(4.0), sq(6.0)], &grid, BurnOptions::default()).unwrap();
        for threads in [1, 0, 8] {
            assert_eq!(burn_par(&[sq(4.0), sq(6.0)], &grid, BurnOptions::default(), threads).unwrap(), want);
        }
        let one = burn(&[sq(4.0)], &grid, BurnOptions::default()).unwrap();
        assert_eq!(burn_par_from(&[sq(4.0)], &grid, BurnOptions::default(), 8, 0).unwrap(), one);
        let none = burn(&[], &grid, BurnOptions::default()).unwrap();
        assert_eq!(burn_par_from(&[], &grid, BurnOptions::default(), 8, 0).unwrap(), none);
        let blobs: [Vec<u8>; 0] = [];
        assert_eq!(burn_wkb_par(&blobs, &grid, BurnOptions::default(), 8).unwrap(), none);
    }

    #[test]
    fn wkb_errors_as_serial() {
        let sq = |x: f64| to_wkb(&Geometry::Polygon(Polygon::new(vec![ring(x, 5.0, 3.0, 12, &[0.5])])), false, 0);
        let blobs = [sq(3.0), sq(6.0)];
        let ok = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
        let bad = GridSpec::new(0.0, 0.0, 10.0, 10.0, 0, 10);
        for (grid, opts) in [
            (bad, BurnOptions::default()),
            (ok, BurnOptions::default().with_id_offset(-1)),
            (ok, BurnOptions::default().with_id_offset(i32::MAX)),
        ] {
            let want = burn_wkb(blobs.iter().map(|b| b.as_slice()), &grid, opts).unwrap_err().to_string();
            assert_eq!(burn_wkb_par_from(&blobs, &grid, opts, 4, 0).unwrap_err().to_string(), want);
        }
    }
}
