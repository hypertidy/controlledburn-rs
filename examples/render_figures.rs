//! Render the book's figures into book/src/img/. Run by hand when a
//! figure changes (`cargo run --release --example render_figures`); the
//! images are committed so the book build never depends on this step.
//! Pure Rust (plotters, bitmap backend), rectangles only: no fonts.

use controlledburn::{
    burn, burn_wkb, materialize, BurnOptions, BurnResult, Coord, EdgePolicy, Geometry, GridSpec, MaterializeOptions,
    PixelFn, Polygon,
};
use plotters::prelude::*;

const OUT: &str = "book/src/img";

fn shade(frac: f64) -> RGBColor {
    // white (0) -> deep orange (1)
    let f = frac.clamp(0.0, 1.0);
    RGBColor((255.0 - 55.0 * f) as u8, (255.0 - 150.0 * f) as u8, (255.0 - 210.0 * f) as u8)
}

/// Paint a BurnResult onto a bitmap, `scale` pixels per cell, shading
/// each cell by its total coverage (sum of fractions over all
/// geometries), computed with the crate's own `materialize`.
fn paint(path: &str, r: &BurnResult, ncol: u32, nrow: u32, scale: u32, grid_lines: bool) {
    let mut cov = vec![f64::NAN; ncol as usize * nrow as usize];
    let opts = MaterializeOptions { fn_: PixelFn::Sum, edge_policy: EdgePolicy::Fraction, threshold: 0.5 };
    materialize(r, &mut cov, ncol, nrow, Some(&vec![1.0; 1 << 16]), &opts).unwrap();

    let (w, h) = (ncol * scale, nrow * scale);
    let root = BitMapBackend::new(path, (w, h)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    for row in 0..nrow {
        for col in 0..ncol {
            let v = cov[(row * ncol + col) as usize];
            if v.is_nan() {
                continue;
            }
            let (x, y) = ((col * scale) as i32, (row * scale) as i32);
            root.draw(&Rectangle::new([(x, y), (x + scale as i32, y + scale as i32)], shade(v).filled())).unwrap();
        }
    }
    if grid_lines {
        let g = RGBColor(180, 180, 180);
        for c in 0..=ncol {
            let x = (c * scale) as i32;
            root.draw(&Rectangle::new([(x, 0), (x + 1, h as i32)], g.filled())).unwrap();
        }
        for rr in 0..=nrow {
            let y = (rr * scale) as i32;
            root.draw(&Rectangle::new([(0, y), (w as i32, y + 1)], g.filled())).unwrap();
        }
    }
    root.present().unwrap();
    println!("wrote {path}");
}

fn main() {
    std::fs::create_dir_all(OUT).unwrap();

    // Chapter 1: the square on the 10 x 10 grid.
    let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
    let ring = vec![
        Coord::new(2.5, 4.5),
        Coord::new(6.5, 4.5),
        Coord::new(6.5, 8.5),
        Coord::new(2.5, 8.5),
        Coord::new(2.5, 4.5),
    ];
    let square = Geometry::Polygon(Polygon::new(vec![ring.clone()]));
    let r = burn(&[square], &grid, BurnOptions::default()).unwrap();
    let path = format!("{OUT}/contract.png");
    // Cells shaded by fraction, grid lines, and the polygon outline.
    {
        let root = BitMapBackend::new(&path, (400, 400)).into_drawing_area();
        root.fill(&WHITE).unwrap();
        let cell = |row: i32, col: i32, c: RGBColor| {
            let (x, y) = ((col - 1) * 40, (row - 1) * 40);
            root.draw(&Rectangle::new([(x, y), (x + 40, y + 40)], c.filled())).unwrap();
        };
        for run in &r.runs {
            for col in run.col_start..=run.col_end {
                cell(run.row, col, shade(1.0));
            }
        }
        for e in &r.edges {
            cell(e.row, e.col, shade(e.fraction as f64));
        }
        let g = RGBColor(180, 180, 180);
        for i in 0..=10 {
            root.draw(&Rectangle::new([(i * 40, 0), (i * 40 + 1, 400)], g.filled())).unwrap();
            root.draw(&Rectangle::new([(0, i * 40), (400, i * 40 + 1)], g.filled())).unwrap();
        }
        let to_px = |c: &Coord| -> (i32, i32) { ((c.x * 40.0) as i32, ((10.0 - c.y) * 40.0) as i32) };
        for s in ring.windows(2) {
            let (x0, y0) = to_px(&s[0]);
            let (x1, y1) = to_px(&s[1]);
            let n = (x1 - x0).abs().max((y1 - y0).abs()).max(1);
            for i in 0..=n {
                let x = x0 + (x1 - x0) * i / n;
                let y = y0 + (y1 - y0) * i / n;
                root.draw(&Rectangle::new([(x - 1, y - 1), (x + 1, y + 1)], BLACK.filled())).unwrap();
            }
        }
        root.present().unwrap();
        println!("wrote {path}");
    }

    // Chapter 2: the world at one degree.
    let bytes = std::fs::read("book/data/ne_110m_countries.wkb").unwrap();
    let n = u32::from_le_bytes(bytes[0..4].try_into().unwrap()) as usize;
    let mut pos = 4;
    let mut blobs = Vec::with_capacity(n);
    for _ in 0..n {
        let len = u32::from_le_bytes(bytes[pos..pos + 4].try_into().unwrap()) as usize;
        pos += 4;
        blobs.push(&bytes[pos..pos + len]);
        pos += len;
    }
    let grid = GridSpec::new(-180.0, -90.0, 180.0, 90.0, 360, 180);
    let r = burn_wkb(blobs.iter().copied(), &grid, BurnOptions::default()).unwrap();
    paint(&format!("{OUT}/world-1deg.png"), &r, 360, 180, 3, false);
}
