# External oracle fixtures (from burnscars)

These files are copied unchanged from
[earth-mover/burnscars](https://github.com/earth-mover/burnscars) at commit
`1d24ef7f14a4dd21396f3d627cdc32796c278d05`, Copyright 2026 The burnscars
authors, Apache License 2.0:

| file | upstream path | what it is |
|---|---|---|
| `cases.csv` | `fixtures/oracle/cases.csv` | 50 cases: name, WKB (hex), grid extent and size |
| `coverage.csv` | `fixtures/oracle/coverage.csv` | per-cell coverage from exactextract 0.3.0 |
| `approx.csv` | `fixtures/oracle/approx.csv` | cells inside by GDAL's rasterize (cell-centre rule), through rasterio |
| `known_diffs.toml` | `oracle/known_diffs.toml` | the cells where burnscars differs from the oracle, each with its reason |
| `corpus.wkt` | `oracle/corpus.wkt` | the hand-written cases |
| `generate.py` | `oracle/generate.py` | the generator: the corpus plus 40 random draws with a fixed seed, kept when valid with area > 0 |

To regenerate, use the upstream `oracle/` project, which pins exactextract,
GDAL and shapely with `uv.lock`: `uv run --locked --project oracle python
oracle/generate.py --out <dir>`.

`tests/oracle.rs` requires cb's difference set for each case to be
exactly the one listed here, less 9 approx cells where cb agrees with
GDAL and burnscars does not (cell centres within about 4e-14 of an edge,
listed in the test). The coverage differences are the 5 cells where
exactextract 0.3.0 gives 0 for a full interior cell; cb and burnscars
both give 1.
