//! Render the book's figures into book/src/img/ as SVG. Run by hand when a
//! figure changes (`cargo run --release --example render_figures`); the
//! images are committed so the book build never depends on this step.
//!
//! Every figure is drawn from a real burn with the crate's own functions;
//! the SVG is written by hand (rectangles, polygons and a little text), so
//! the renderer needs no dependencies.

use std::fmt::Write as _;

use controlledburn::{
    burn, burn_wkb, crop, materialize, parse_wkb, BurnOptions, BurnResult, Coord, EdgePolicy, Geometry, GridSpec,
    MaterializeOptions, PixelFn, Polygon,
};

const OUT: &str = "book/src/img";

// Coverage: one hue, white (0) to deep orange (1), as in the book so far.
const COVER: (u8, u8, u8) = (200, 105, 45);
// Line length: one hue, white to blue.
const LENGTH: (u8, u8, u8) = (42, 120, 214);
// Diverging, approx minus coverage: blue (approx counts more), red (less).
const MORE: (u8, u8, u8) = (42, 120, 214);
const LESS: (u8, u8, u8) = (227, 73, 72);
const INK: &str = "#333333";
const MUTED: &str = "#777777";
const GRID: &str = "#d4d4d4";
const FONT: &str = "font-family=\"system-ui, -apple-system, 'Segoe UI', Roboto, sans-serif\"";

fn mix(c: (u8, u8, u8), f: f64) -> String {
    // 32 steps: invisible, and same-coloured neighbours merge into one path.
    let f = (f.clamp(0.0, 1.0) * 32.0).round() / 32.0;
    let ch = |v: u8| (255.0 - (255.0 - v as f64) * f).round() as u8;
    format!("#{:02x}{:02x}{:02x}", ch(c.0), ch(c.1), ch(c.2))
}

struct Svg {
    w: f64,
    h: f64,
    body: String,
    clips: usize,
}

impl Svg {
    fn new(w: f64, h: f64) -> Self {
        let mut s = Svg { w, h, body: String::new(), clips: 0 };
        s.rect(0.0, 0.0, w, h, "#ffffff");
        s
    }

    fn rect(&mut self, x: f64, y: f64, w: f64, h: f64, fill: &str) {
        let _ = writeln!(self.body, r#"<rect x="{x:.2}" y="{y:.2}" width="{w:.2}" height="{h:.2}" fill="{fill}"/>"#);
    }

    fn outline(&mut self, [x, y, w, h]: [f64; 4], stroke: &str, width: f64, dash: bool) {
        let d = if dash { r#" stroke-dasharray="6 4""# } else { "" };
        let _ = writeln!(
            self.body,
            r#"<rect x="{x:.2}" y="{y:.2}" width="{w:.2}" height="{h:.2}" fill="none" stroke="{stroke}" stroke-width="{width}"{d}/>"#
        );
    }

    fn line(&mut self, x0: f64, y0: f64, x1: f64, y1: f64, stroke: &str, width: f64) {
        let _ = writeln!(
            self.body,
            r#"<line x1="{x0:.2}" y1="{y0:.2}" x2="{x1:.2}" y2="{y1:.2}" stroke="{stroke}" stroke-width="{width}"/>"#
        );
    }

    fn quad(&mut self, pts: &[(f64, f64)], fill: &str, stroke: &str) {
        let p: Vec<String> = pts.iter().map(|(x, y)| format!("{x:.2},{y:.2}")).collect();
        let _ = writeln!(
            self.body,
            r#"<polygon points="{}" fill="{fill}" stroke="{stroke}" stroke-width="0.6"/>"#,
            p.join(" ")
        );
    }

    /// Rings as one even-odd path, clipped to a rectangle.
    fn rings(&mut self, rings: &[Vec<(f64, f64)>], fill: &str, stroke: &str, width: f64, clip: [f64; 4]) {
        let mut d = String::new();
        for r in rings {
            for (i, (x, y)) in r.iter().enumerate() {
                let _ = write!(d, "{}{x:.1},{y:.1}", if i == 0 { "M" } else { "L" });
            }
            d.push('Z');
        }
        if d.is_empty() {
            return;
        }
        self.clips += 1;
        let id = format!("c{}", self.clips);
        let [x, y, w, h] = clip;
        let _ = writeln!(
            self.body,
            r#"<clipPath id="{id}"><rect x="{x:.2}" y="{y:.2}" width="{w:.2}" height="{h:.2}"/></clipPath>"#
        );
        let _ = writeln!(
            self.body,
            r#"<path d="{d}" fill="{fill}" fill-rule="evenodd" stroke="{stroke}" stroke-width="{width}" stroke-linejoin="round" clip-path="url(#{id})"/>"#
        );
    }

    fn dot(&mut self, x: f64, y: f64, r: f64, fill: &str) {
        let _ = writeln!(
            self.body,
            r##"<circle cx="{x:.2}" cy="{y:.2}" r="{r}" fill="{fill}" stroke="#ffffff" stroke-width="1.5"/>"##
        );
    }

    fn text(&mut self, x: f64, y: f64, s: &str, size: f64, anchor: &str, fill: &str) {
        let _ = writeln!(
            self.body,
            r#"<text x="{x:.2}" y="{y:.2}" font-size="{size}" text-anchor="{anchor}" fill="{fill}" {FONT}>{s}</text>"#
        );
    }

    fn save(&self, name: &str) {
        let path = format!("{OUT}/{name}");
        let svg = format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}\" height=\"{h}\" viewBox=\"0 0 {w} {h}\">\n{}</svg>\n",
            self.body,
            w = self.w,
            h = self.h
        );
        std::fs::write(&path, svg).unwrap();
        println!("wrote {path} ({} KB)", std::fs::metadata(&path).unwrap().len() / 1024);
    }
}

/// A grid drawn at `(x0, y0)` with `s` pixels a cell.
#[derive(Clone, Copy)]
struct Panel {
    x0: f64,
    y0: f64,
    s: f64,
    grid: GridSpec,
}

impl Panel {
    fn w(&self) -> f64 {
        self.grid.ncol as f64 * self.s
    }
    fn h(&self) -> f64 {
        self.grid.nrow as f64 * self.s
    }
    fn cell(&self, row: i32, col: i32) -> (f64, f64) {
        (self.x0 + col as f64 * self.s, self.y0 + row as f64 * self.s)
    }
    fn xy(&self, c: &Coord) -> (f64, f64) {
        let g = &self.grid;
        (self.x0 + (c.x - g.xmin) / g.dx() * self.s, self.y0 + (g.ymax - c.y) / g.dy() * self.s)
    }
    fn clip(&self) -> [f64; 4] {
        [self.x0, self.y0, self.w(), self.h()]
    }

    fn paint_cell(&self, svg: &mut Svg, row: i32, col: i32, fill: &str) {
        let (x, y) = self.cell(row, col);
        svg.rect(x, y, self.s, self.s, fill);
    }

    /// Cell spans `(row, col_start, col_end, colour)`, drawn as one path
    /// per colour with neighbours of the same colour merged.
    fn spans(&self, svg: &mut Svg, mut spans: Vec<(i32, i32, i32, String)>) {
        spans.sort_by_key(|a| (a.0, a.1));
        let mut merged: Vec<(i32, i32, i32, String)> = Vec::new();
        for sp in spans {
            match merged.last_mut() {
                Some(m) if m.0 == sp.0 && m.2 == sp.1 && m.3 == sp.3 => m.2 = sp.2,
                _ => merged.push(sp),
            }
        }
        let mut by: std::collections::BTreeMap<String, String> = Default::default();
        for (row, c0, c1, colour) in merged {
            let (x, y) = self.cell(row, c0);
            let d = by.entry(colour).or_default();
            let _ = write!(
                d,
                "M{},{}h{}v{}h-{}z",
                num(x),
                num(y),
                num((c1 - c0) as f64 * self.s),
                num(self.s),
                num((c1 - c0) as f64 * self.s)
            );
        }
        for (colour, d) in by {
            let _ = writeln!(svg.body, r#"<path fill="{colour}" d="{d}"/>"#);
        }
    }

    /// Runs solid, edges shaded by fraction, in one hue.
    fn records(&self, svg: &mut Svg, r: &BurnResult, hue: (u8, u8, u8)) {
        self.records_at(svg, r, hue, 1.0);
    }

    /// As `records`, at a fraction of full strength.
    fn records_at(&self, svg: &mut Svg, r: &BurnResult, hue: (u8, u8, u8), k: f64) {
        let mut sp: Vec<_> = r.runs.iter().map(|run| (run.row, run.col_start, run.col_end, mix(hue, k))).collect();
        sp.extend(r.edges.iter().map(|e| (e.row, e.col, e.col + 1, mix(hue, k * e.fraction as f64))));
        self.spans(svg, sp);
    }

    /// A dense buffer, NaN left blank.
    fn pixels(&self, svg: &mut Svg, px: &[f64], fill: impl Fn(f64) -> Option<String>) {
        let nc = self.grid.ncol as usize;
        let sp = px
            .iter()
            .enumerate()
            .filter(|(_, v)| !v.is_nan())
            .filter_map(|(i, &v)| fill(v).map(|f| ((i / nc) as i32, (i % nc) as i32, (i % nc) as i32 + 1, f)))
            .collect();
        self.spans(svg, sp);
    }

    /// True if a geometry's bounding box meets the grid's.
    fn meets(&self, g: &Geometry) -> bool {
        let polys: Vec<&Polygon> = match g {
            Geometry::Polygon(p) => vec![p],
            Geometry::MultiPolygon(ps) => ps.iter().collect(),
            _ => return false,
        };
        let gr = &self.grid;
        polys
            .iter()
            .flat_map(|p| p.rings.iter().flatten())
            .any(|c| c.x >= gr.xmin - 5.0 && c.x <= gr.xmax + 5.0 && c.y >= gr.ymin - 5.0 && c.y <= gr.ymax + 5.0)
    }

    fn lines(&self, svg: &mut Svg, colour: &str, width: f64) {
        for c in 0..=self.grid.ncol {
            let x = self.x0 + c as f64 * self.s;
            svg.line(x, self.y0, x, self.y0 + self.h(), colour, width);
        }
        for r in 0..=self.grid.nrow {
            let y = self.y0 + r as f64 * self.s;
            svg.line(self.x0, y, self.x0 + self.w(), y, colour, width);
        }
    }

    fn frame(&self, svg: &mut Svg) {
        svg.outline([self.x0, self.y0, self.w(), self.h()], MUTED, 1.0, false);
    }

    fn geometry(&self, svg: &mut Svg, g: &Geometry, fill: &str, stroke: &str, width: f64) {
        let polys: Vec<&Polygon> = match g {
            Geometry::Polygon(p) => vec![p],
            Geometry::MultiPolygon(ps) => ps.iter().collect(),
            _ => vec![],
        };
        for p in polys {
            let rings: Vec<Vec<(f64, f64)>> = p.rings.iter().map(|r| r.iter().map(|c| self.xy(c)).collect()).collect();
            svg.rings(&rings, fill, stroke, width, self.clip());
        }
    }

    fn title(&self, svg: &mut Svg, s: &str) {
        svg.text(self.x0, self.y0 - 10.0, s, 14.0, "start", INK);
    }
}

/// A coordinate with at most one decimal, no trailing zero.
fn num(v: f64) -> String {
    let r = (v * 10.0).round() / 10.0;
    if r.fract() == 0.0 {
        format!("{}", r as i64)
    } else {
        format!("{r:.1}")
    }
}

fn square() -> Vec<Coord> {
    vec![Coord::new(2.5, 4.5), Coord::new(6.5, 4.5), Coord::new(6.5, 8.5), Coord::new(2.5, 8.5), Coord::new(2.5, 4.5)]
}

fn countries() -> Vec<Vec<u8>> {
    let bytes = std::fs::read("book/data/ne_110m_countries.wkb").unwrap();
    let n = u32::from_le_bytes(bytes[0..4].try_into().unwrap()) as usize;
    let mut pos = 4;
    let mut blobs = Vec::with_capacity(n);
    for _ in 0..n {
        let len = u32::from_le_bytes(bytes[pos..pos + 4].try_into().unwrap()) as usize;
        pos += 4;
        blobs.push(bytes[pos..pos + len].to_vec());
        pos += len;
    }
    blobs
}

fn dense(r: &BurnResult, grid: &GridSpec, opts: &MaterializeOptions) -> Vec<f64> {
    let mut px = vec![f64::NAN; grid.ncol as usize * grid.nrow as usize];
    let ones = vec![1.0; 1 << 16];
    materialize(r, &mut px, grid.ncol, grid.nrow, Some(&ones), opts).unwrap();
    px
}

fn fmt_frac(f: f32) -> String {
    let s = format!("{f:.2}");
    s.trim_start_matches('0').to_string()
}

/// The contract: the square on the 10 x 10 grid, runs solid, edges shaded
/// and labelled with their fraction.
fn contract() {
    let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
    let ring = square();
    let g = Geometry::Polygon(Polygon::new(vec![ring]));
    let r = burn(std::slice::from_ref(&g), &grid, BurnOptions::default()).unwrap();
    let p = Panel { x0: 30.0, y0: 30.0, s: 40.0, grid };
    let mut svg = Svg::new(p.w() + 50.0, p.h() + 60.0);
    p.records(&mut svg, &r, COVER);
    p.lines(&mut svg, GRID, 1.0);
    p.geometry(&mut svg, &g, "none", INK, 2.0);
    // Each label nudged off the outline, towards the square's middle
    // (row 3, col 4).
    for e in &r.edges {
        let (x, y) = p.cell(e.row, e.col);
        let dx = (4 - e.col).signum() as f64 * 0.22 * p.s;
        let dy = (3 - e.row).signum() as f64 * 0.22 * p.s;
        svg.text(x + p.s / 2.0 + dx, y + p.s / 2.0 + 4.0 + dy, &fmt_frac(e.fraction), 10.0, "middle", INK);
    }
    for i in 0..10 {
        let c = p.x0 + (i as f64 + 0.5) * p.s;
        svg.text(c, p.y0 - 8.0, &i.to_string(), 11.0, "middle", MUTED);
        svg.text(p.x0 - 8.0, p.y0 + (i as f64 + 0.5) * p.s + 4.0, &i.to_string(), 11.0, "end", MUTED);
    }
    svg.text(p.x0 + p.w() / 2.0, p.y0 + p.h() + 22.0, "col", 11.0, "middle", MUTED);
    svg.save("contract.svg");
}

/// Coverage and approx on a coastline, and where they differ.
fn coverage_vs_approx(blobs: &[Vec<u8>]) {
    // Britain and Ireland at a quarter degree.
    let grid = GridSpec::new(-11.0, 49.5, 2.5, 59.5, 54, 40);
    let wkb = || blobs.iter().map(|b| b.as_slice());
    let cov = burn_wkb(wkb(), &grid, BurnOptions::coverage()).unwrap();
    let app = burn_wkb(wkb(), &grid, BurnOptions::approx()).unwrap();
    let sum = MaterializeOptions::new(PixelFn::Sum, EdgePolicy::Fraction);
    let pc = dense(&cov, &grid, &sum);
    let pa = dense(&app, &grid, &sum);
    let geoms: Vec<Geometry> = blobs.iter().map(|b| parse_wkb(b).unwrap()).collect();

    let s = 6.0;
    let gap = 24.0;
    let w = grid.ncol as f64 * s;
    let panels: Vec<Panel> = (0..3).map(|i| Panel { x0: 10.0 + i as f64 * (w + gap), y0: 30.0, s, grid }).collect();
    let mut svg = Svg::new(20.0 + 3.0 * w + 2.0 * gap, 40.0 + grid.nrow as f64 * s);

    panels[0].records(&mut svg, &cov, COVER);
    panels[0].title(&mut svg, "coverage: exact fractions");
    panels[1].records(&mut svg, &app, COVER);
    panels[1].title(&mut svg, "approx: cell centre inside");
    let diff: Vec<f64> = pa
        .iter()
        .zip(&pc)
        .map(|(a, c)| {
            let v = |x: &f64| if x.is_nan() { 0.0 } else { x.min(1.0) };
            v(a) - v(c)
        })
        .collect();
    panels[2].pixels(&mut svg, &diff, |d| {
        if d.abs() < 1e-9 {
            None
        } else if d > 0.0 {
            Some(mix(MORE, d.max(0.15)))
        } else {
            Some(mix(LESS, (-d).max(0.15)))
        }
    });
    panels[2].title(&mut svg, "approx minus coverage");
    for p in &panels {
        for g in geoms.iter().filter(|g| p.meets(g)) {
            p.geometry(&mut svg, g, "none", INK, 0.7);
        }
        p.frame(&mut svg);
    }
    svg.save("coverage-vs-approx.svg");
}

/// Lines: length in each cell. Points: the cell of each point.
fn lines_and_points() {
    let grid = GridSpec::new(0.0, 0.0, 12.0, 8.0, 12, 8);
    let line = Geometry::LineString(vec![
        Coord::new(0.6, 1.2),
        Coord::new(3.7, 6.4),
        Coord::new(7.2, 2.3),
        Coord::new(11.4, 7.1),
    ]);
    let pts = vec![
        Coord::new(1.5, 6.5),
        Coord::new(3.2, 2.7),
        Coord::new(3.8, 2.2),
        Coord::new(6.4, 6.6),
        Coord::new(8.7, 1.4),
        Coord::new(10.2, 4.8),
    ];
    let points = Geometry::MultiPoint(pts.clone());
    let rl = burn(std::slice::from_ref(&line), &grid, BurnOptions::default()).unwrap();
    let rp = burn(std::slice::from_ref(&points), &grid, BurnOptions::default()).unwrap();

    let s = 30.0;
    let gap = 30.0;
    let a = Panel { x0: 10.0, y0: 30.0, s, grid };
    let b = Panel { x0: 10.0 + a.w() + gap, ..a };
    let mut svg = Svg::new(20.0 + 2.0 * a.w() + gap, 40.0 + a.h());

    let diag = 2f64.sqrt();
    for l in &rl.lines {
        a.paint_cell(&mut svg, l.row, l.col, &mix(LENGTH, l.length as f64 / diag));
    }
    a.lines(&mut svg, GRID, 1.0);
    if let Geometry::LineString(c) = &line {
        for w in c.windows(2) {
            let (x0, y0) = a.xy(&w[0]);
            let (x1, y1) = a.xy(&w[1]);
            svg.line(x0, y0, x1, y1, INK, 2.0);
        }
    }
    a.frame(&mut svg);
    a.title(&mut svg, &format!("a line: {} records, length in each cell", rl.lines.len()));

    for p in &rp.points {
        b.paint_cell(&mut svg, p.row, p.col, &mix(COVER, 0.55));
    }
    b.lines(&mut svg, GRID, 1.0);
    for c in &pts {
        let (x, y) = b.xy(c);
        svg.dot(x, y, 4.5, INK);
    }
    b.frame(&mut svg);
    b.title(&mut svg, &format!("{} points: {} records, one per point", pts.len(), rp.points.len()));
    svg.save("lines-and-points.svg");
}

/// The world at one degree, cells shaded by coverage.
fn world(blobs: &[Vec<u8>]) {
    let grid = GridSpec::new(-180.0, -90.0, 180.0, 90.0, 360, 180);
    let r = burn_wkb(blobs.iter().map(|b| b.as_slice()), &grid, BurnOptions::default()).unwrap();
    let px = dense(&r, &grid, &MaterializeOptions::new(PixelFn::Sum, EdgePolicy::Fraction));
    let p = Panel { x0: 0.0, y0: 0.0, s: 2.0, grid };
    let mut svg = Svg::new(p.w(), p.h());
    p.pixels(&mut svg, &px, |v| Some(mix(COVER, v)));
    svg.save("world-1deg.svg");
}

/// The world at one degree, burned once; one 60 x 45 degree tile cropped
/// and materialized, the rest never made dense.
fn materialize_tile(blobs: &[Vec<u8>]) {
    let grid = GridSpec::new(-180.0, -90.0, 180.0, 90.0, 360, 180);
    let r = burn_wkb(blobs.iter().map(|b| b.as_slice()), &grid, BurnOptions::default()).unwrap();
    let p = Panel { x0: 10.0, y0: 10.0, s: 2.0, grid };
    let mut svg = Svg::new(p.w() + 20.0, p.h() + 20.0);
    // The tables, faint: each cell's records summed, as materialize would.
    let mut sum: std::collections::HashMap<(i32, i32), f64> = Default::default();
    for e in &r.edges {
        *sum.entry((e.row, e.col)).or_default() += e.fraction as f64;
    }
    let mut sp: Vec<_> = r.runs.iter().map(|run| (run.row, run.col_start, run.col_end, mix(COVER, 0.22))).collect();
    sp.extend(sum.into_iter().map(|((row, col), f)| (row, col, col + 1, mix(COVER, 0.22 * f.min(1.0)))));
    p.spans(&mut svg, sp);
    // The tile.
    let target = [0.0, 0.0, 60.0, 45.0];
    let (tile, sub) = crop(&r, &grid, target).unwrap();
    let px = dense(&tile, &sub, &MaterializeOptions::new(PixelFn::Sum, EdgePolicy::Fraction));
    let tp =
        Panel { x0: p.xy(&Coord::new(sub.xmin, 0.0)).0, y0: p.xy(&Coord::new(0.0, sub.ymax)).1, s: p.s, grid: sub };
    svg.rect(tp.x0, tp.y0, tp.w(), tp.h(), "#ffffff");
    tp.pixels(&mut svg, &px, |v| Some(mix(COVER, v)));
    // Tile boundaries.
    for i in 0..=6 {
        let x = p.x0 + i as f64 * 60.0 * p.s;
        svg.line(x, p.y0, x, p.y0 + p.h(), GRID, 1.0);
    }
    for j in 0..=4 {
        let y = p.y0 + j as f64 * 45.0 * p.s;
        svg.line(p.x0, y, p.x0 + p.w(), y, GRID, 1.0);
    }
    svg.outline([tp.x0, tp.y0, tp.w(), tp.h()], INK, 2.0, false);
    svg.save("materialize-tile.svg");
}

/// Crop: a window cut from the burn, its indices re-based to 0.
fn crop_figure() {
    let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
    let g = Geometry::Polygon(Polygon::new(vec![square()]));
    let r = burn(std::slice::from_ref(&g), &grid, BurnOptions::default()).unwrap();
    let (tile, sub) = crop(&r, &grid, [3.0, 5.0, 7.0, 9.0]).unwrap();

    let s = 30.0;
    let a = Panel { x0: 40.0, y0: 50.0, s, grid };
    let b = Panel { x0: a.x0 + a.w() + 90.0, y0: a.y0 + 3.0 * s, s, grid: sub };
    let mut svg = Svg::new(b.x0 + b.w() + 90.0, a.y0 + a.h() + 20.0);

    a.records(&mut svg, &r, COVER);
    a.lines(&mut svg, GRID, 1.0);
    let (wx, wy) = a.xy(&Coord::new(3.0, 9.0));
    svg.outline([wx, wy, 4.0 * s, 4.0 * s], INK, 2.5, true);
    a.title(&mut svg, "burn, with the window [3, 5, 7, 9]");
    b.records(&mut svg, &tile, COVER);
    b.lines(&mut svg, GRID, 1.0);
    b.frame(&mut svg);
    b.title(&mut svg, "crop: 4 x 4, from (0, 0)");
    for i in 0..10 {
        svg.text(a.x0 + (i as f64 + 0.5) * s, a.y0 + a.h() + 14.0, &i.to_string(), 10.0, "middle", MUTED);
        svg.text(a.x0 - 6.0, a.y0 + (i as f64 + 0.5) * s + 4.0, &i.to_string(), 10.0, "end", MUTED);
    }
    for i in 0..4 {
        svg.text(b.x0 + (i as f64 + 0.5) * s, b.y0 + b.h() + 14.0, &i.to_string(), 10.0, "middle", MUTED);
        svg.text(b.x0 - 6.0, b.y0 + (i as f64 + 0.5) * s + 4.0, &i.to_string(), 10.0, "end", MUTED);
    }
    svg.line(wx + 4.0 * s + 8.0, wy + 2.0 * s, b.x0 - 22.0, b.y0 + 2.0 * s, MUTED, 1.0);
    svg.save("crop.svg");
}

/// merge_parts: two parts meeting mid-cell, off and on.
fn merge_parts() {
    let grid = GridSpec::new(0.0, 0.0, 9.0, 4.0, 9, 4);
    let rect = |x0: f64, x1: f64| {
        Polygon::new(vec![vec![
            Coord::new(x0, 0.6),
            Coord::new(x1, 0.6),
            Coord::new(x1, 3.4),
            Coord::new(x0, 3.4),
            Coord::new(x0, 0.6),
        ]])
    };
    let g = Geometry::MultiPolygon(vec![rect(1.3, 4.5), rect(4.5, 7.7)]);
    let off = burn(std::slice::from_ref(&g), &grid, BurnOptions::default().with_merge_parts(false)).unwrap();
    let on = burn(std::slice::from_ref(&g), &grid, BurnOptions::default()).unwrap();

    let s = 36.0;
    let a = Panel { x0: 10.0, y0: 30.0, s, grid };
    let b = Panel { y0: a.y0 + a.h() + 50.0, ..a };
    let mut svg = Svg::new(a.w() + 20.0, b.y0 + b.h() + 10.0);
    for (p, r, name) in [(&a, &off, "merge_parts off"), (&b, &on, "merge_parts on (default)")] {
        // Every record of a cell side by side, so a cell with two records
        // shows two.
        let mut cells: std::collections::BTreeMap<(i32, i32), Vec<f64>> = Default::default();
        for run in &r.runs {
            for c in run.col_start..run.col_end {
                cells.entry((run.row, c)).or_default().push(1.0);
            }
        }
        for e in &r.edges {
            cells.entry((e.row, e.col)).or_default().push(e.fraction as f64);
        }
        for (&(row, col), fs) in &cells {
            let (x, y) = p.cell(row, col);
            let w = p.s / fs.len() as f64;
            for (i, f) in fs.iter().enumerate() {
                svg.rect(x + i as f64 * w, y, w, p.s, &mix(COVER, *f));
                if i > 0 {
                    svg.line(x + i as f64 * w, y, x + i as f64 * w, y + p.s, "#ffffff", 2.0);
                }
            }
            if fs.len() > 1 {
                svg.text(x + p.s * 0.8, y + p.s * 0.35, &fs.len().to_string(), 10.0, "middle", INK);
            }
        }
        p.lines(&mut svg, GRID, 1.0);
        p.geometry(&mut svg, &g, "none", INK, 1.5);
        let n = r.runs.len() + r.edges.len();
        p.title(&mut svg, &format!("{name}: {n} records"));
    }
    svg.save("merge-parts.svg");
}

/// A rotated grid: burn in index space, read the records in world space.
fn grids() {
    let (ncol, nrow) = (10u32, 7u32);
    let theta = 25f64.to_radians();
    let (ox, oy) = (2.0, 0.0);
    // world = origin + col * u + row * v; u along a row, v down a column.
    let u = (theta.cos(), theta.sin());
    let v = (theta.sin(), -theta.cos());
    let world = |col: f64, row: f64| Coord::new(ox + col * u.0 + row * v.0, oy + col * u.1 + row * v.1);
    let index = |c: &Coord| {
        let (dx, dy) = (c.x - ox, c.y - oy);
        // u and v are orthonormal.
        (dx * u.0 + dy * u.1, dx * v.0 + dy * v.1)
    };
    let poly = [world(1.4, 1.1), world(6.2, 0.6), world(8.7, 3.2), world(6.9, 6.3), world(2.3, 5.4), world(1.4, 1.1)];
    // Index space: row negated, on unit cells.
    let lattice = GridSpec::new(0.0, -(nrow as f64), ncol as f64, 0.0, ncol, nrow);
    let ring: Vec<Coord> = poly
        .iter()
        .map(|c| {
            let (col, row) = index(c);
            Coord::new(col, -row)
        })
        .collect();
    let gi = Geometry::Polygon(Polygon::new(vec![ring]));
    let r = burn(std::slice::from_ref(&gi), &lattice, BurnOptions::default()).unwrap();

    let s = 32.0;
    let a = Panel { x0: 20.0, y0: 40.0, s, grid: lattice };
    // World panel: bounding box of the grid's corners.
    let corners = [world(0.0, 0.0), world(ncol as f64, 0.0), world(0.0, nrow as f64), world(ncol as f64, nrow as f64)];
    let (wx0, wx1) = corners.iter().fold((f64::MAX, f64::MIN), |(a, b), c| (a.min(c.x), b.max(c.x)));
    let (wy0, wy1) = corners.iter().fold((f64::MAX, f64::MIN), |(a, b), c| (a.min(c.y), b.max(c.y)));
    let bx = a.x0 + a.w() + 60.0;
    let to = |c: &Coord| (bx + (c.x - wx0) * s, a.y0 + (wy1 - c.y) * s);
    let mut svg = Svg::new(bx + (wx1 - wx0) * s + 20.0, a.y0 + ((wy1 - wy0) * s).max(a.h()) + 20.0);

    a.records(&mut svg, &r, COVER);
    a.lines(&mut svg, GRID, 1.0);
    a.geometry(&mut svg, &gi, "none", INK, 1.5);
    a.title(&mut svg, "index space: burn on unit cells");

    let mut frac = vec![0.0f64; (ncol * nrow) as usize];
    for run in &r.runs {
        for c in run.col_start..run.col_end {
            frac[(run.row as u32 * ncol + c as u32) as usize] = 1.0;
        }
    }
    for e in &r.edges {
        frac[(e.row as u32 * ncol + e.col as u32) as usize] = e.fraction as f64;
    }
    for row in 0..nrow {
        for col in 0..ncol {
            let (c, rr) = (col as f64, row as f64);
            let q: Vec<(f64, f64)> = [world(c, rr), world(c + 1.0, rr), world(c + 1.0, rr + 1.0), world(c, rr + 1.0)]
                .iter()
                .map(to)
                .collect();
            let f = frac[(row * ncol + col) as usize];
            svg.quad(&q, &if f > 0.0 { mix(COVER, f) } else { "#ffffff".to_string() }, GRID);
        }
    }
    let pw: Vec<(f64, f64)> = poly.iter().map(to).collect();
    svg.rings(&[pw], "none", INK, 1.5, [0.0, 0.0, svg.w, svg.h]);
    let (tx, ty) = to(&world(0.0, 0.0));
    svg.dot(tx, ty, 4.0, INK);
    svg.text(tx - 8.0, ty - 6.0, "(0, 0)", 11.0, "end", MUTED);
    let (lx, ly) = a.cell(0, 0);
    svg.dot(lx, ly, 4.0, INK);
    svg.text(bx, a.y0 - 10.0, "world: the same records on the rotated grid", 14.0, "start", INK);
    svg.save("grids.svg");
}

fn main() {
    std::fs::create_dir_all(OUT).unwrap();
    let blobs = countries();
    contract();
    coverage_vs_approx(&blobs);
    lines_and_points();
    world(&blobs);
    materialize_tile(&blobs);
    crop_figure();
    merge_parts();
    grids();
}
