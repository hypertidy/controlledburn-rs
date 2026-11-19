//! Sparse rasterization output: the four-table contract.
//!
//! Polygon -> `runs` (interior RLE) + `edges` (boundary coverage fractions).
//! Line    -> `lines` (length in cell, CRS units).
//! Point   -> `points` (no measure column; implicit 1).
//!
//! Row and column indices are 0-based and row 0 is the top row; a run's
//! `col_end` is exclusive, so `col_end - col_start` is its length and
//! `buf[col_start..col_end]` is its slice. Schemas are
//! type-pure: each table's measure column (or absence thereof) means
//! exactly one thing.

/// A single interior run: fully covered cells `col_start..col_end`
/// (half-open) on `row`.
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GridRun {
    pub row: i32,
    pub col_start: i32,
    pub col_end: i32,
    pub id: i32,
}

/// A single polygon-boundary cell. `fraction` is the dimensionless
/// coverage fraction in (0, 1): area(polygon intersect cell) / area(cell).
///
/// `fraction` is computed in `f64` and stored as `f32`, good to about 6e-8
/// of the cell; this keeps the record at 16 bytes. For `f64` fractions over
/// a window use [`CoverageWindow`](crate::CoverageWindow).
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GridEdge {
    pub row: i32,
    pub col: i32,
    pub fraction: f32,
    pub id: i32,
}

/// A single line cell. `length` is the absolute length of the line within
/// the cell, in CRS units.
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GridLine {
    pub row: i32,
    pub col: i32,
    pub length: f32,
    pub id: i32,
}

/// A single point cell.
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GridPoint {
    pub row: i32,
    pub col: i32,
    pub id: i32,
}

/// A non-fatal problem encountered for one input geometry.
/// `geom_index` is the 0-based position in the input (the same as `id`).
#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Note {
    pub geom_index: i32,
    pub message: String,
}

/// The result of a burn.
///
/// # Reading the tables
///
/// - A cell that is in no record is 0 for every id: nothing of any
///   geometry falls in it.
/// - For one id a cell has at most one record: a run or an edge (or a
///   line, or points). Across ids a cell can have one record for each, so
///   two polygons that share a border both have an edge in the cells it
///   crosses, with fractions that add up to 1 where they tile the cell.
/// - Each table is grouped by `id`, in input order. Within an id, `runs`,
///   `edges` and `lines` are in row-major order: by `row`, then by column
///   (`col_start` for a run). `points` are in the order the geometry
///   lists them, one record per point, so a cell repeats when two points
///   fall in it.
/// - With [`BurnOptions::merge_parts`](crate::BurnOptions::merge_parts)
///   off (and so with `BurnOptions::parity`), a multi-part geometry has
///   records for each part in turn, each part in row-major order, and a
///   cell two parts cross has a record from each.
///
/// The order holds for every burn function, and [`crop`](crate::crop),
/// [`slice_tiles`](crate::slice_tiles) and
/// [`slice_tiles_par`](crate::slice_tiles_par) keep it. Results put
/// together afterwards (concatenated, or folded by a front end) get it
/// back from [`merge_cells`](crate::merge_cells).
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BurnResult {
    pub runs: Vec<GridRun>,
    pub edges: Vec<GridEdge>,
    pub lines: Vec<GridLine>,
    pub points: Vec<GridPoint>,
    pub notes: Vec<Note>,
}

impl BurnResult {
    /// Total number of cells covered by `runs` (sum of run lengths).
    pub fn run_cells(&self) -> u64 {
        self.runs.iter().map(|r| (r.col_end - r.col_start) as u64).sum()
    }

    /// Sum of run cells plus edge fractions: the covered area in cell units.
    pub fn covered_cells(&self) -> f64 {
        self.run_cells() as f64 + self.edges.iter().map(|e| e.fraction as f64).sum::<f64>()
    }

    /// Total line length across all line cells, in CRS units.
    pub fn line_length(&self) -> f64 {
        self.lines.iter().map(|l| l.length as f64).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.runs.is_empty() && self.edges.is_empty() && self.lines.is_empty() && self.points.is_empty()
    }

    pub(crate) fn note(&mut self, id: i32, message: impl Into<String>) {
        self.notes.push(Note { geom_index: id, message: message.into() });
    }
}
