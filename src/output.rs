//! Sparse rasterization output: the four-table contract.
//!
//! Polygon -> `runs` (interior RLE) + `edges` (boundary coverage fractions).
//! Line    -> `lines` (length in cell, CRS units).
//! Point   -> `points` (no measure column; implicit 1).
//!
//! All row/col indices are 1-based and row 1 is the top row. Schemas are
//! type-pure: each table's measure column (or absence thereof) means
//! exactly one thing.

/// A single interior run: fully covered cells `[col_start, col_end]` on `row`.
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
/// `geom_index` is the 1-based position in the input.
#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Note {
    pub geom_index: i32,
    pub message: String,
}

/// The result of a burn.
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
        self.runs
            .iter()
            .map(|r| (r.col_end - r.col_start + 1) as u64)
            .sum()
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

    pub(crate) fn note(&mut self, geom_index: usize, message: impl Into<String>) {
        self.notes.push(Note {
            geom_index: geom_index as i32,
            message: message.into(),
        });
    }
}
