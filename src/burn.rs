//! Public entry points: `burn` and `burn_wkb`.

use std::borrow::Cow;

use crate::approx::process_polygon_approx;
use crate::coalesce::coalesce_sorted;
use crate::error::BurnError;
use crate::geometry::{Geometry, Polygon};
use crate::grid::{Grid, GridSpec};
use crate::line::process_line;
use crate::output::BurnResult;
use crate::point::process_point;
use crate::polygon::process_polygon;
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

/// Options for a burn.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BurnOptions {
    pub mode: BurnMode,
    /// Merge adjacent runs on a row (default `true`). With `false` the
    /// `runs` table is record-for-record identical to the C++ core, which
    /// emits one run per fully covered boundary cell.
    pub coalesce_runs: bool,
}

impl Default for BurnOptions {
    fn default() -> Self {
        BurnOptions { mode: BurnMode::Coverage, coalesce_runs: true }
    }
}

impl BurnOptions {
    /// Coverage mode with coalescing on.
    pub fn coverage() -> Self {
        Self::default()
    }

    /// Approx mode with coalescing on.
    pub fn approx() -> Self {
        BurnOptions { mode: BurnMode::Approx, coalesce_runs: true }
    }

    /// Record-for-record parity with the C++ core: coalescing off.
    pub fn parity(mode: BurnMode) -> Self {
        BurnOptions { mode, coalesce_runs: false }
    }

    pub fn with_mode(mut self, mode: BurnMode) -> Self {
        self.mode = mode;
        self
    }

    pub fn with_coalesce_runs(mut self, on: bool) -> Self {
        self.coalesce_runs = on;
        self
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

fn process_geometry(g: &Geometry, full_grid: &Grid<0>, gs: &GridSpec, geom_id: i32, out: &mut BurnResult, opts: &BurnOptions) -> Result<(), String> {
    let run_start = out.runs.len();
    for poly in g.polygons() {
        match opts.mode {
            BurnMode::Approx => process_polygon_approx(poly, gs, geom_id, &mut out.runs),
            BurnMode::Coverage => process_polygon(poly, full_grid, geom_id, &mut out.runs, &mut out.edges)?,
        }
        if opts.coalesce_runs && out.runs.len() > run_start {
            let mut tail = out.runs.split_off(run_start);
            coalesce_sorted(&mut tail);
            out.runs.extend(tail);
        }
    }
    for line in g.lines() {
        process_line(line, full_grid, geom_id, &mut out.lines)?;
    }
    for pt in g.points() {
        process_point(pt, full_grid, geom_id, &mut out.points);
    }
    Ok(())
}

fn burn_one(g: &Geometry, k: usize, full_grid: &Grid<0>, gs: &GridSpec, out: &mut BurnResult, opts: &BurnOptions) {
    if g.is_empty() {
        return;
    }
    if !g.is_finite() {
        out.note(k + 1, "skipped geometry with non-finite coordinates");
        return;
    }
    let g = normalize(g);
    if let Err(e) = process_geometry(&g, full_grid, gs, (k + 1) as i32, out, opts) {
        out.note(k + 1, format!("error processing geometry: {e}"));
    }
}

/// Burn a set of geometries onto the grid. Geometry `k` (0-based) is
/// assigned id `k + 1` in the output tables.
///
/// Empty geometries are skipped silently; geometries with non-finite
/// coordinates and geometries the engine cannot process are skipped with
/// a note. Unclosed polygon rings are closed. The only error is an
/// invalid `GridSpec`.
pub fn burn(geoms: &[Geometry], grid: &GridSpec, opts: impl Into<BurnOptions>) -> Result<BurnResult, BurnError> {
    grid.validate()?;
    let opts = opts.into();
    let full_grid = grid.bounded();
    let mut out = BurnResult::default();
    for (k, g) in geoms.iter().enumerate() {
        burn_one(g, k, &full_grid, grid, &mut out, &opts);
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
    let full_grid = grid.bounded();
    let mut out = BurnResult::default();
    for (k, blob) in wkb.into_iter().enumerate() {
        if blob.is_empty() {
            continue;
        }
        match parse_wkb(blob) {
            Ok(g) => burn_one(&g, k, &full_grid, grid, &mut out, &opts),
            Err(e) => out.note(k + 1, format!("failed to parse WKB: {e}")),
        }
    }
    Ok(out)
}
