// Derived from exactextract side.h
// Copyright (c) 2018 ISciences, LLC. Apache License 2.0.

/// A side of an axis-aligned box.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Side {
    Left,
    Right,
    Top,
    Bottom,
}
