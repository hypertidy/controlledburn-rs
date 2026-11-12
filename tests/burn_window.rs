//! `burn_window` against the whole-grid burn: for random polygons (with
//! holes, multi-part), lattice-aligned rectangles, lines along cell walls
//! and points on cell corners, on a grid whose cell size is not exact in
//! binary (at the origin and far from it), and tile shapes that usually do not divide it, every window's
//! burn equals the whole burn's tile from `slice_tiles`, record for
//! record, in both modes and with and without coalescing and merging.

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

/// The whole burn's tile, back on the full grid's indices.
fn shifted(t: &Tile) -> BurnResult {
    let (dr, dc) = (t.row_off as i32, t.col_off as i32);
    let mut r = t.result.clone();
    r.runs.iter_mut().for_each(|x| {
        x.row += dr;
        x.col_start += dc;
        x.col_end += dc;
    });
    r.edges.iter_mut().for_each(|x| {
        x.row += dr;
        x.col += dc;
    });
    r.lines.iter_mut().for_each(|x| {
        x.row += dr;
        x.col += dc;
    });
    r.points.iter_mut().for_each(|x| {
        x.row += dr;
        x.col += dc;
    });
    r
}

/// 17 x 11 units in 13 x 7 cells, at an origin away from zero, as a
/// projected grid's is.
fn grid(o: (f64, f64)) -> GridSpec {
    GridSpec::new(o.0, o.1, o.0 + 17.0, o.1 + 11.0, NCOL, NROW)
}

#[allow(clippy::too_many_arguments)]
fn geoms(o: (f64, f64), cx: f64, cy: f64, r: f64, k: usize, wobble: &[f64], c: [u32; 4]) -> Vec<Geometry> {
    let grid = grid(o);
    let (cx, cy) = (cx + o.0, cy + o.1);
    let ring = |cx: f64, cy: f64, r, k, w| ring(cx, cy, r, k, w);
    let (dx, dy) = (grid.dx(), grid.dy());
    let x = |c: u32| grid.xmin + c as f64 * dx;
    let y = |r: u32| grid.ymax - r as f64 * dy;
    let (c0, c1) = (c[0].min(c[1]), c[0].max(c[1]) + 1);
    let (r0, r1) = (c[2].min(c[3]), c[2].max(c[3]) + 1);
    let rect = vec![
        Coord::new(x(c0), y(r1)),
        Coord::new(x(c1), y(r1)),
        Coord::new(x(c1), y(r0)),
        Coord::new(x(c0), y(r0)),
        Coord::new(x(c0), y(r1)),
    ];
    let mut hole = ring(cx, cy, r * 0.3, k, wobble);
    hole.reverse();
    vec![
        Geometry::Polygon(Polygon::new(vec![ring(cx, cy, r, k, wobble)])),
        Geometry::MultiPolygon(vec![
            Polygon::new(vec![ring(o.0 + (cy - o.1), o.1 + (cx - o.0) * 0.6, r * 0.8, k, wobble), hole]),
            Polygon::new(vec![ring(o.0 + (cx - o.0) * 0.5, o.1 + (cy - o.1) * 0.5, r * 0.4, k, wobble)]),
        ]),
        // corners on the lattice: walls of the window and of the grid
        Geometry::Polygon(Polygon::new(vec![rect.clone()])),
        Geometry::LineString(ring(o.0 + (cy - o.1), o.1 + (cx - o.0), r * 0.7, k, wobble)),
        // along walls and through corners
        Geometry::MultiLineString(vec![
            rect,
            vec![Coord::new(x(0), y(c[2])), Coord::new(x(NCOL), y(c[2]))],
            vec![Coord::new(x(c[0]), y(0)), Coord::new(x(c[0]), y(NROW))],
        ]),
        // on corners and walls, so on window seams and the grid's edges
        Geometry::MultiPoint(
            (0..=NCOL).map(|c| Coord::new(x(c), y(c % (NROW + 1)))).chain([Coord::new(x(c0), cy)]).collect(),
        ),
    ]
}

proptest! {
    #[test]
    fn windows_equal_the_whole_burns_tiles(
        cx in 0.0..17.0f64, cy in 0.0..11.0f64, r in 0.5..9.0f64,
        k in 3usize..12, wobble in prop::collection::vec(0.0..1.0f64, 12),
        c in prop::array::uniform4(0u32..7), tc in 1u32..16, tr in 1u32..9,
        approx in any::<bool>(), coalesce in any::<bool>(), merge in any::<bool>(),
        o in prop_oneof![Just((0.0, 0.0)), (-1e6..1e6f64, -1e7..1e7f64)],
    ) {
        let grid = grid(o);
        let c = [c[0] * 2 % (NCOL + 1), c[1] * 2 % (NCOL + 1), c[2], c[3]];
        let c = [c[0].min(NCOL - 1), c[1].min(NCOL - 1), c[2].min(NROW - 1), c[3].min(NROW - 1)];
        let g = geoms(o, cx, cy, r, k, &wobble, c);
        let mode = if approx { BurnMode::Approx } else { BurnMode::Coverage };
        let opts = BurnOptions::default().with_mode(mode).with_coalesce_runs(coalesce).with_merge_parts(merge);
        let whole = burn(&g, &grid, opts).unwrap();
        let tiles = slice_tiles(&whole, &grid, tc, tr);
        for ti in 0..NROW.div_ceil(tr) {
            for tj in 0..NCOL.div_ceil(tc) {
                let rows = ti * tr..((ti + 1) * tr).min(NROW);
                let cols = tj * tc..((tj + 1) * tc).min(NCOL);
                let mut w = burn_window(&g, &grid, rows, cols, opts).unwrap();
                prop_assert_eq!(&w.notes, &whole.notes);
                w.notes.clear();
                let want = tiles
                    .iter()
                    .find(|t| t.tile_row == ti && t.tile_col == tj)
                    .map_or_else(BurnResult::default, shifted);

                prop_assert_eq!(&w, &want, "tile ({}, {}) of {} x {}", ti, tj, tc, tr);
            }
        }
    }
}

#[test]
fn a_window_must_be_in_the_grid() {
    let grid = grid((0.0, 0.0));
    let opts = BurnOptions::default();
    assert!(burn_window(&[], &grid, 0..NROW, 0..NCOL, opts).is_ok());
    for (rows, cols) in [(0..0, 0..1), (0..1, 3..3), (0..NROW + 1, 0..1), (0..1, 2..NCOL + 1)] {
        assert!(matches!(burn_window(&[], &grid, rows, cols, opts), Err(BurnError::InvalidGrid(_))));
    }
}

/// A vertex at lattice position (a, b) / 4 of a cell: on a node, a wall
/// or inside a cell, as `a` and `b` fall.
fn lattice(grid: &GridSpec, a: i32, b: i32) -> Coord {
    Coord::new(grid.xmin + a as f64 / 4.0 * grid.dx(), grid.ymax - b as f64 / 4.0 * grid.dy())
}

proptest! {
    /// Long lines on a larger grid with small windows, so most of each
    /// line is trimmed: vertices on nodes and walls, segments along walls
    /// (staircases), and lines that leave the grid and come back.
    #[test]
    fn line_windows_equal_the_whole_burns_tiles(
        steps in prop::collection::vec((-2i32..3, -2i32..3, 0u8..4), 2..60),
        start in (0i32..188, 0i32..124), tc in 1u32..9, tr in 1u32..9, merge in any::<bool>(),
        o in prop_oneof![Just((0.0, 0.0)), (-1e6..1e6f64, -1e7..1e7f64)],
    ) {
        let (ncol, nrow) = (47, 31);
        let grid = GridSpec::new(o.0, o.1, o.0 + 47.0 * 0.3, o.1 + 31.0 * 0.7, ncol, nrow);
        let (mut a, mut b) = (start.0 - 4, start.1 - 4);
        let mut pts = vec![lattice(&grid, a, b)];
        for &(da, db, kind) in &steps {
            // 0: to a node, 1: along a wall, 2: anywhere on the quarter
            // lattice, 3: a long jump
            match kind {
                0 => { a = (a / 4 + da) * 4; b = (b / 4 + db) * 4; }
                1 => { if da.abs() > db.abs() { a = (a / 4 + da) * 4; b = b / 4 * 4 } else { b = (b / 4 + db) * 4; a = a / 4 * 4 } }
                2 => { a += da * 3; b += db * 5; }
                _ => { a += da * 37; b += db * 23; }
            }
            pts.push(lattice(&grid, a, b));
        }
        let half = pts.len() / 2;
        let g = vec![
            Geometry::LineString(pts.clone()),
            Geometry::MultiLineString(vec![pts[..half.max(2)].to_vec(), pts[half.saturating_sub(1)..].to_vec()]),
        ];
        let opts = BurnOptions::default().with_merge_parts(merge);
        let whole = burn(&g, &grid, opts).unwrap();
        let tiles = slice_tiles(&whole, &grid, tc, tr);
        for ti in 0..nrow.div_ceil(tr) {
            for tj in 0..ncol.div_ceil(tc) {
                let rows = ti * tr..((ti + 1) * tr).min(nrow);
                let cols = tj * tc..((tj + 1) * tc).min(ncol);
                let mut w = burn_window(&g, &grid, rows, cols, opts).unwrap();
                w.notes.clear();
                let want = tiles
                    .iter()
                    .find(|t| t.tile_row == ti && t.tile_col == tj)
                    .map_or_else(BurnResult::default, shifted);
                prop_assert_eq!(&w, &want, "tile ({}, {}) of {} x {}", ti, tj, tc, tr);
            }
        }
    }
}
