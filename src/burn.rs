//! Public entry points: `burn` and `burn_wkb`.

use std::borrow::Cow;

use crate::approx::process_polygon_approx;
use crate::coalesce::coalesce_sorted_from;
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

fn process_geometry(
    g: &Geometry,
    full_grid: &Grid<0>,
    gs: &GridSpec,
    geom_id: i32,
    out: &mut BurnResult,
    opts: &BurnOptions,
    sc: &mut CoverageScratch,
) -> Result<(), String> {
    let run_start = out.runs.len();
    let edge_start = out.edges.len();
    let line_start = out.lines.len();
    let polys = g.polygons();
    note_short_rings(polys, geom_id, out);
    for poly in polys {
        match opts.mode {
            BurnMode::Approx => process_polygon_approx(poly, gs, geom_id, &mut out.runs),
            BurnMode::Coverage => process_polygon(poly, full_grid, geom_id, &mut out.runs, &mut out.edges, sc)?,
        }
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
        process_line(line, full_grid, geom_id, &mut out.lines)?;
    }
    if opts.merge_parts && g.lines().len() > 1 {
        merge_part_lines(out, line_start);
    }
    for pt in g.points() {
        process_point(pt, full_grid, geom_id, &mut out.points);
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

fn burn_one(
    g: &Geometry,
    id: i32,
    full_grid: &Grid<0>,
    gs: &GridSpec,
    out: &mut BurnResult,
    opts: &BurnOptions,
    sc: &mut CoverageScratch,
) {
    if g.is_empty() {
        return;
    }
    if !g.is_finite() {
        out.note(id, "skipped geometry with non-finite coordinates");
        return;
    }
    let g = normalize(g);
    if let Err(e) = process_geometry(&g, full_grid, gs, id, out, opts, sc) {
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
        burn_one(g, opts.id(k)?, &full_grid, grid, &mut out, &opts, &mut sc);
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
            Ok(g) => burn_one(&g, id, &full_grid, grid, &mut out, &opts, &mut sc),
            Err(e) => out.note(id, format!("failed to parse WKB: {e}")),
        }
    }
    Ok(out)
}
