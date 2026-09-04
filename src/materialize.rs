//! Optional consumer: burn a sparse `BurnResult` into a caller-provided
//! pixel buffer with per-pixel reduction functions (fasterize semantics).
//!
//! The buffer is row-major, row 0 (top) first, `ncol * nrow` values,
//! caller-initialised (conventionally to NaN as the background). Cells
//! never touched are left untouched.
//!
//! Semantic note: fasterize burns a cell when the cell CENTRE is inside
//! the polygon; a Coverage-mode result carries exact fractions instead.
//! Interior runs agree; `edge_policy` decides boundary cells:
//! `Fraction` writes `value * fraction` (area conserving), `Threshold`
//! includes the cell iff `fraction >= threshold` (0.5 approximates, but
//! does not equal, the centre rule; use `BurnMode::Approx` for that).

use crate::error::BurnError;
use crate::output::BurnResult;

/// Reduction applied when several geometries touch one pixel.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum PixelFn {
    /// Keep the first value written.
    First,
    /// Overwrite with each new value.
    #[default]
    Last,
    Sum,
    Min,
    Max,
    /// Number of geometries touching the pixel (ignores value).
    Count,
    /// 1 if touched.
    Any,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum EdgePolicy {
    /// value * coverage fraction, combined per `PixelFn`.
    Fraction,
    /// All-or-nothing at `threshold`.
    #[default]
    Threshold,
}

#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MaterializeOptions {
    pub fn_: PixelFn,
    pub edge_policy: EdgePolicy,
    pub threshold: f64,
}

impl Default for MaterializeOptions {
    fn default() -> Self {
        MaterializeOptions { fn_: PixelFn::Last, edge_policy: EdgePolicy::Threshold, threshold: 0.5 }
    }
}

#[inline]
fn apply(px: &mut f64, value: f64, f: PixelFn) {
    let empty = px.is_nan();
    match f {
        PixelFn::First => {
            if empty {
                *px = value;
            }
        }
        PixelFn::Last => *px = value,
        PixelFn::Sum => *px = if empty { value } else { *px + value },
        PixelFn::Min => {
            if empty || value < *px {
                *px = value;
            }
        }
        PixelFn::Max => {
            if empty || value > *px {
                *px = value;
            }
        }
        PixelFn::Count => *px = if empty { 1.0 } else { *px + 1.0 },
        PixelFn::Any => *px = 1.0,
    }
}

/// Materialize polygon output (runs + edges) into `buffer`.
///
/// `values` maps geometry id to burn value: the value for id `k` is
/// `values[k]`. Pass `None` to burn the id itself.
pub fn materialize(
    result: &BurnResult,
    buffer: &mut [f64],
    ncol: u32,
    nrow: u32,
    values: Option<&[f64]>,
    opts: &MaterializeOptions,
) -> Result<(), BurnError> {
    if ncol == 0 || nrow == 0 {
        return Err(BurnError::InvalidBuffer("invalid dimensions".into()));
    }
    let expected = ncol as usize * nrow as usize;
    if buffer.len() != expected {
        return Err(BurnError::InvalidBuffer(format!(
            "buffer length {} does not match ncol * nrow = {expected}",
            buffer.len()
        )));
    }
    let ncol_i = ncol as i32;
    let nrow_i = nrow as i32;

    let value_of = |id: i32| -> Result<f64, BurnError> {
        match values {
            None => Ok(id as f64),
            Some(v) => v.get(id as usize).copied().ok_or(BurnError::IdOutOfRange { id, values: v.len() }),
        }
    };
    let idx = |row: i32, col: i32| -> usize { row as usize * ncol as usize + col as usize };

    for r in &result.runs {
        if r.row < 0 || r.row >= nrow_i {
            continue;
        }
        let c0 = r.col_start.max(0);
        let c1 = r.col_end.min(ncol_i);
        let v = value_of(r.id)?;
        for c in c0..c1 {
            apply(&mut buffer[idx(r.row, c)], v, opts.fn_);
        }
    }

    for e in &result.edges {
        if e.row < 0 || e.row >= nrow_i || e.col < 0 || e.col >= ncol_i {
            continue;
        }
        let v = value_of(e.id)?;
        match opts.edge_policy {
            EdgePolicy::Threshold => {
                if e.fraction as f64 >= opts.threshold {
                    apply(&mut buffer[idx(e.row, e.col)], v, opts.fn_);
                }
            }
            EdgePolicy::Fraction => apply(&mut buffer[idx(e.row, e.col)], v * e.fraction as f64, opts.fn_),
        }
    }
    Ok(())
}
