#!/bin/sh
# Build the crate into a dedicated target dir (so only one rlib matches)
# and run every code block in the book against it.
set -e
cd "$(dirname "$0")/.."
cargo build --lib --all-features --target-dir target/book
D=target/book/debug
if [ -d "$D/deps" ]; then
  # Cargo before Rust 1.100: every rlib in debug/deps.
  LIBS="$D/deps"
else
  # Cargo from Rust 1.100: the crate's rlib is copied to debug/, and each
  # dependency's (rlib, or shared library for a proc macro) is in
  # debug/build/<crate>/<hash>/out. Leave out the crate's own copy there,
  # or rustdoc finds two candidates.
  LIBS="$D"
  for f in $(find "$D/build" \( -name '*.rlib' -o -name '*.so' -o -name '*.dylib' -o -name '*.dll' \) ! -name 'libcontrolledburn-*'); do
    LIBS="$LIBS,$(dirname "$f")"
  done
fi
CONTROLLEDBURN_BOOK_DATA="$PWD/book/data" mdbook test book -L "$LIBS"
