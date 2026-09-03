#!/bin/sh
# Build the crate into a dedicated target dir (so only one rlib matches)
# and run every code block in the book against it.
set -e
cd "$(dirname "$0")/.."
cargo build --lib --all-features --target-dir target/book
CONTROLLEDBURN_BOOK_DATA="$PWD/book/data" mdbook test book -L target/book/debug/deps
