//! `crop` against slicing the full dense raster: for random polygons,
//! lines and points, and random windows (cell-aligned or not, inside,
//! straddling or beyond the grid), materializing the cropped tables on
//! the returned sub-grid equals the matching cells of the full raster.

use controlledburn::*;
use proptest::prelude::*;

fn dense(r: &BurnResult, ncol: u32, nrow: u32) -> Vec<f64> {
    let mut buf = vec![0.0; (ncol * nrow) as usize];
    let opts = MaterializeOptions { fn_: PixelFn::Sum, edge_policy: EdgePolicy::Fraction, threshold: 0.5 };
    materialize(r, &mut buf, ncol, nrow, None, &opts).unwrap();
    buf
}

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

proptest! {
    #[test]
    fn crop_equals_slice(
        cx in 0.0..17.0f64, cy in 0.0..11.0f64, r in 0.5..9.0f64,
        k in 3usize..12, wobble in prop::collection::vec(0.0..1.0f64, 12),
        w in prop::array::uniform4(-4.0..21.0f64),
    ) {
        // a non-dyadic grid: 17 x 11 units in 13 x 7 cells
        let grid = GridSpec::new(0.0, 0.0, 17.0, 11.0, 13, 7);
        let geoms = vec![
            Geometry::Polygon(Polygon::new(vec![ring(cx, cy, r, k, &wobble)])),
            Geometry::LineString(ring(cy, cx, r * 0.7, k, &wobble)),
            Geometry::MultiPoint(ring(cx, cy, r * 0.5, k, &wobble)),
        ];
        let full = burn(&geoms, &grid, BurnOptions::default()).unwrap();
        let fd = dense(&full, grid.ncol, grid.nrow);
        // lines and points too: sum their measures on top of coverage
        let target = [w[0].min(w[1]), w[2].min(w[3]) * 11.0 / 21.0, w[0].max(w[1]), w[2].max(w[3]) * 11.0 / 21.0];
        match crop(&full, &grid, target) {
            None => {
                // no cell overlaps: the window misses the grid or is thinner
                // than the snapping tolerance
                let (dx, dy) = (grid.dx(), grid.dy());
                let c0 = ((target[0] / dx) + 1e-8).floor().max(0.0);
                let c1 = ((target[2] / dx) - 1e-8).ceil().min(13.0);
                let r0 = (((11.0 - target[3]) / dy) + 1e-8).floor().max(0.0);
                let r1 = (((11.0 - target[1]) / dy) - 1e-8).ceil().min(7.0);
                prop_assert!(c0 >= c1 || r0 >= r1);
            }
            Some((sub, g)) => {
                let (dx, dy) = (grid.dx(), grid.dy());
                let c0 = ((g.xmin - grid.xmin) / dx).round() as usize;
                let r0 = ((grid.ymax - g.ymax) / dy).round() as usize;
                prop_assert!((g.dx() - dx).abs() < 1e-9 && (g.dy() - dy).abs() < 1e-9);
                prop_assert!(c0 + g.ncol as usize <= 13 && r0 + g.nrow as usize <= 7);
                for run in &sub.runs {
                    prop_assert!(run.col_start >= 0 && run.col_start < run.col_end && run.col_end <= g.ncol as i32);
                    prop_assert!(run.row >= 0 && run.row < g.nrow as i32);
                }
                let sd = dense(&sub, g.ncol, g.nrow);
                for rr in 0..g.nrow as usize {
                    for cc in 0..g.ncol as usize {
                        let a = sd[rr * g.ncol as usize + cc];
                        let b = fd[(r0 + rr) * 13 + c0 + cc];
                        prop_assert_eq!(a.to_bits(), b.to_bits(), "cell ({}, {})", rr, cc);
                    }
                }
                // and nothing inside the window was lost
                let cells: f64 = sd.iter().sum();
                let mut want = 0.0;
                for rr in 0..g.nrow as usize {
                    for cc in 0..g.ncol as usize {
                        want += fd[(r0 + rr) * 13 + c0 + cc];
                    }
                }
                prop_assert_eq!(cells.to_bits(), want.to_bits());
            }
        }
    }
}
