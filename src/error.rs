//! Error types.

use std::fmt;

/// The only fatal error a burn can raise: an invalid grid or buffer.
/// Everything else is reported per geometry in `BurnResult::notes`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BurnError {
    /// `ncol` or `nrow` is zero, or the extent is empty or inverted.
    InvalidGrid(String),
    /// `materialize` was given a buffer of the wrong length or invalid dimensions.
    InvalidBuffer(String),
    /// `materialize` saw a geometry id beyond the supplied values.
    IdOutOfRange { id: i32, values: usize },
}

impl fmt::Display for BurnError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BurnError::InvalidGrid(m) => write!(f, "invalid grid: {m}"),
            BurnError::InvalidBuffer(m) => write!(f, "materialize: {m}"),
            BurnError::IdOutOfRange { id, values } => {
                write!(f, "materialize: geometry id {id} exceeds values size {values}")
            }
        }
    }
}

impl std::error::Error for BurnError {}

/// A WKB parse failure. `Display` renders the same message text as the
/// C++ core so notes stay comparable.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum WkbError {
    TooShort,
    BadByteOrder(u8),
    Truncated { at: usize },
    MismatchedPart { found: u32 },
    GeometryCollection,
    Unsupported { type_code: u32 },
}

impl fmt::Display for WkbError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            WkbError::TooShort => write!(f, "WKB too short"),
            WkbError::BadByteOrder(_) => write!(f, "WKB invalid byte order flag"),
            WkbError::Truncated { at } => write!(f, "WKB truncated at byte {at}"),
            WkbError::MismatchedPart { found } => {
                write!(f, "WKB multi-geometry contains mismatched part type {found}")
            }
            WkbError::GeometryCollection => write!(
                f,
                "GeometryCollection is not supported (mixed dimensions break weight semantics); split into homogeneous groups"
            ),
            WkbError::Unsupported { type_code } => write!(
                f,
                "WKB unsupported geometry type {type_code} (curved types must be linearised upstream)"
            ),
        }
    }
}

impl std::error::Error for WkbError {}
