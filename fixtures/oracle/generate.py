"""Oracle fixtures for burnscars.

Coverage oracle: exactextract (cell_id, coverage). Approx oracle: GDAL rasterize through
rasterio (centre rule). GEOS (shapely) is used only as a black box: is_valid, area, WKB.
"""

import argparse
import math
from pathlib import Path

import numpy as np
import rasterio.features
import rasterio.transform
import shapely
from exactextract import exact_extract
from exactextract.raster import NumPyRasterSource


def read_corpus(path):
    cases = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        name, wkt, grid = line.split("|")
        xmin, ymin, xmax, ymax, ncol, nrow = grid.split()
        cases.append((name, shapely.from_wkt(wkt), (float(xmin), float(ymin), float(xmax), float(ymax), int(ncol), int(nrow))))
    return cases


def random_cases(seed=20261014, n=40):
    rng = np.random.default_rng(seed)
    out = []
    for i in range(n):
        ncol, nrow = int(rng.integers(5, 60)), int(rng.integers(5, 60))
        dx, dy = 10.0 / ncol, 10.0 / nrow
        k = int(rng.integers(3, 30))
        cx, cy = rng.uniform(2, 8, 2)
        ang = np.sort(rng.uniform(0, 2 * math.pi, k))
        rad = rng.uniform(0.3, 4.0, k)
        xs, ys = cx + rad * np.cos(ang), cy + rad * np.sin(ang)
        if i % 3 == 0:  # vertices on the lattice
            xs, ys = np.round(xs / dx) * dx, np.round(ys / dy) * dy
        if i % 3 == 1:  # vertices a few ulp from the lattice
            xs = np.nextafter(np.round(xs / dx) * dx, np.inf)
        g = shapely.Polygon(np.column_stack([xs, ys]))
        if g.is_valid and g.area > 0:
            out.append((f"random_{i}", g, (0.0, 0.0, 10.0, 10.0, ncol, nrow)))
    return out


# exactextract 0.3.0 gives 0 for some interior cells when vertices are on or 1 ulp off the lattice
# (a per-cell check with GEOS area shows it). oracle/known_diffs.toml lists these cases.
AREA_CHECK_FAILS = {"random_1", "random_4", "random_28", "random_39"}


def coverage(name, g, grid):
    xmin, ymin, xmax, ymax, ncol, nrow = grid
    src = NumPyRasterSource(np.zeros((nrow, ncol)), xmin, ymin, xmax, ymax)
    feats = [{"type": "Feature", "properties": {}, "geometry": shapely.geometry.mapping(g)}]
    res = exact_extract(src, feats, ["cell_id", "coverage"], output="geojson")
    props = res[0]["properties"]
    cells = [(int(c) // ncol, int(c) % ncol, float(f)) for c, f in zip(props["cell_id"], props["coverage"])]
    # Self-check of the cell_id meaning: the sum of coverage is the clipped area.
    cell_area = (xmax - xmin) / ncol * (ymax - ymin) / nrow
    want = g.intersection(shapely.box(xmin, ymin, xmax, ymax)).area
    got = sum(f for _, _, f in cells) * cell_area
    ok = abs(got - want) <= 1e-6 * max(want, 1.0)
    # An exempt case must still fail, so that an oracle fix shows here.
    assert ok != (name in AREA_CHECK_FAILS), (name, got, want)
    return cells


def approx(g, grid):
    xmin, ymin, xmax, ymax, ncol, nrow = grid
    t = rasterio.transform.from_bounds(xmin, ymin, xmax, ymax, ncol, nrow)
    a = rasterio.features.rasterize([(g, 1)], out_shape=(nrow, ncol), transform=t, all_touched=False, fill=0, dtype="uint8")
    return [(int(r), int(c)) for r, c in zip(*np.nonzero(a))]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True)
    p.add_argument("--corpus", default=str(Path(__file__).with_name("corpus.wkt")))
    p.add_argument("--no-random", action="store_true")
    a = p.parse_args()
    cases = read_corpus(a.corpus) + ([] if a.no_random else random_cases())
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "cases.csv", "w") as fc, open(out / "coverage.csv", "w") as fv, open(out / "approx.csv", "w") as fa:
        fc.write("case,wkb_hex,xmin,ymin,xmax,ymax,ncol,nrow\n")
        fv.write("case,row,col,fraction\n")
        fa.write("case,row,col\n")
        for name, g, grid in cases:
            assert g.is_valid, name
            fc.write(f"{name},{shapely.to_wkb(g, hex=True)},{','.join(repr(v) for v in grid)}\n")
            for r, c, f in coverage(name, g, grid):
                if f > 0:
                    fv.write(f"{name},{r},{c},{f:.9g}\n")
            for r, c in approx(g, grid):
                fa.write(f"{name},{r},{c}\n")


if __name__ == "__main__":
    main()
