//! `slice_tiles` against `crop` and against the whole dense raster: for
//! random polygons (with holes), lines and points, on a non-dyadic grid
//! with tile shapes that usually do not divide it, every tile equals
//! `crop` to its window record for record, the empty tiles are exactly
//! the ones left out, and the materialized tiles put back together equal
//! the materialized whole grid bit for bit.

use controlledburn::*;
use proptest::prelude::*;

const NCOL: u32 = 13;
const NROW: u32 = 7;

fn ring(cx: f64, cy: f64, r: f64, k: usize, wobble: &[f64]) -> Vec<Coord> {
    let mut v: Vec<Coord> = (0..k)
        .map(|i| {
            let a = i as f64 / k as f64 * std::f64::consts::TAU;
            let rr = r * (0.6 + 0.4 * wobble[i % wobble.len()]);
            Coord::new(cx + rr * a.cos(), cy + rr * a.sin())
        })
        .collect();
    v.push(v[0]);
    v
}

fn dense(r: &BurnResult, ncol: u32, nrow: u32, fn_: PixelFn, edge_policy: EdgePolicy) -> Vec<f64> {
    let mut buf = vec![f64::NAN; (ncol * nrow) as usize];
    let opts = MaterializeOptions { fn_, edge_policy, threshold: 0.5 };
    materialize(r, &mut buf, ncol, nrow, None, &opts).unwrap();
    buf
}

fn tables_eq(a: &BurnResult, b: &BurnResult) -> bool {
    a.runs == b.runs && a.edges == b.edges && a.lines == b.lines && a.points == b.points
}

proptest! {
    #[test]
    fn tiles_equal_crop_and_whole(
        cx in 0.0..17.0f64, cy in 0.0..11.0f64, r in 0.5..9.0f64,
        k in 3usize..12, wobble in prop::collection::vec(0.0..1.0f64, 12),
        tc in 1u32..16, tr in 1u32..9, approx in any::<bool>(),
    ) {
        // 17 x 11 units in 13 x 7 cells; the hole keeps the multi-part
        // and coalescing paths in play
        let grid = GridSpec::new(0.0, 0.0, 17.0, 11.0, NCOL, NROW);
        let shell = ring(cx, cy, r, k, &wobble);
        let mut hole = ring(cx, cy, r * 0.3, k, &wobble);
        hole.reverse();
        let geoms = vec![
            Geometry::Polygon(Polygon::new(vec![shell])),
            Geometry::MultiPolygon(vec![
                Polygon::new(vec![ring(cy, cx * 0.6, r * 0.8, k, &wobble), hole]),
                Polygon::new(vec![ring(cx * 0.5, cy * 0.5, r * 0.4, k, &wobble)]),
            ]),
            Geometry::LineString(ring(cy, cx, r * 0.7, k, &wobble)),
            // points on cell corners and edges, so on tile seams too
            Geometry::MultiPoint((0..=NCOL).map(|c| Coord::new(c as f64 * 17.0 / 13.0, (c % 8) as f64 * 11.0 / 7.0)).collect()),
        ];
        let opts = if approx { BurnOptions::approx() } else { BurnOptions::default() };
        let full = burn(&geoms, &grid, opts).unwrap();
        let tiles = slice_tiles(&full, &grid, tc, tr);

        // every tile against crop, and the empty tiles are the ones left out
        let (ntc, ntr) = (NCOL.div_ceil(tc), NROW.div_ceil(tr));
        let mut seen = 0usize;
        for ti in 0..ntr {
            for tj in 0..ntc {
                let (r0, c0) = (ti * tr, tj * tc);
                let (r1, c1) = ((r0 + tr).min(NROW), (c0 + tc).min(NCOL));
                let ext = [grid.xmin + c0 as f64 * grid.dx(), grid.ymax - r1 as f64 * grid.dy(),
                           grid.xmin + c1 as f64 * grid.dx(), grid.ymax - r0 as f64 * grid.dy()];
                let (cropped, cg) = crop(&full, &grid, ext).unwrap();
                prop_assert_eq!((cg.ncol, cg.nrow), (c1 - c0, r1 - r0));
                let empty = cropped.runs.is_empty() && cropped.edges.is_empty()
                    && cropped.lines.is_empty() && cropped.points.is_empty();
                let t = tiles.iter().find(|t| t.tile_row == ti && t.tile_col == tj);
                match t {
                    None => prop_assert!(empty, "tile ({}, {}) missing", ti, tj),
                    Some(t) => {
                        seen += 1;
                        prop_assert!(!empty);
                        prop_assert_eq!((t.row_off, t.col_off), (r0, c0));
                        prop_assert_eq!((t.grid.ncol, t.grid.nrow), (c1 - c0, r1 - r0));
                        prop_assert!(tables_eq(&t.result, &cropped), "tile ({}, {}) differs from crop", ti, tj);
                    }
                }
            }
        }
        prop_assert_eq!(seen, tiles.len());
        // row-major order
        prop_assert!(tiles.windows(2).all(|w| (w[0].tile_row, w[0].tile_col) < (w[1].tile_row, w[1].tile_col)));
        // the last tiles end exactly on the grid's outer edges
        for t in &tiles {
            if t.col_off + t.grid.ncol == NCOL { prop_assert_eq!(t.grid.xmax, grid.xmax); }
            if t.row_off + t.grid.nrow == NROW { prop_assert_eq!(t.grid.ymin, grid.ymin); }
        }

        // put the tiles back together and compare with the whole, bit for bit
        for (fn_, pol) in [(PixelFn::Sum, EdgePolicy::Fraction), (PixelFn::Last, EdgePolicy::Fraction),
                           (PixelFn::Count, EdgePolicy::Threshold), (PixelFn::Max, EdgePolicy::Fraction)] {
            let whole = dense(&full, NCOL, NROW, fn_, pol);
            let mut assembled = vec![f64::NAN; (NCOL * NROW) as usize];
            for t in &tiles {
                let d = dense(&t.result, t.grid.ncol, t.grid.nrow, fn_, pol);
                for rr in 0..t.grid.nrow {
                    for cc in 0..t.grid.ncol {
                        assembled[((t.row_off + rr) * NCOL + t.col_off + cc) as usize] = d[(rr * t.grid.ncol + cc) as usize];
                    }
                }
            }
            for (i, (a, b)) in assembled.iter().zip(&whole).enumerate() {
                prop_assert_eq!(a.to_bits(), b.to_bits(), "{:?} cell {}", fn_, i);
            }
        }
    }
}

#[test]
fn one_tile_is_the_whole_result() {
    let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
    let sq = common_square();
    let full = burn(&[sq], &grid, BurnOptions::default()).unwrap();
    let tiles = slice_tiles(&full, &grid, 10, 10);
    assert_eq!(tiles.len(), 1);
    assert_eq!(tiles[0].grid, grid);
    assert!(tables_eq(&tiles[0].result, &full));
    // a tile larger than the grid is the same single tile
    let big = slice_tiles(&full, &grid, 64, 64);
    assert_eq!(big.len(), 1);
    assert_eq!(big[0].grid, grid);
}

#[test]
fn empty_burn_has_no_tiles() {
    let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
    assert!(slice_tiles(&BurnResult::default(), &grid, 3, 3).is_empty());
}

#[test]
#[should_panic(expected = "tile shape")]
fn zero_tile_shape_panics() {
    let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
    slice_tiles(&BurnResult::default(), &grid, 0, 3);
}

fn common_square() -> Geometry {
    Geometry::Polygon(Polygon::new(vec![vec![
        Coord::new(2.5, 2.5),
        Coord::new(7.5, 2.5),
        Coord::new(7.5, 7.5),
        Coord::new(2.5, 7.5),
        Coord::new(2.5, 2.5),
    ]]))
}
