#![allow(dead_code)]

use controlledburn::{BurnResult, Coord, Geometry, GridEdge, GridLine, GridPoint, GridRun, GridSpec, Polygon};
use serde_json::Value;

pub fn fixtures_dir() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures")
}

pub fn load_golden() -> Vec<Value> {
    let text = std::fs::read_to_string(fixtures_dir().join("controlledburn-golden.json")).unwrap();
    serde_json::from_str(&text).unwrap()
}

fn coord(v: &Value) -> Coord {
    Coord::new(v[0].as_f64().unwrap(), v[1].as_f64().unwrap())
}

fn seq(v: &Value) -> Vec<Coord> {
    v.as_array().unwrap().iter().map(coord).collect()
}

/// Geometry from the golden JSON's dumped C++ struct
/// (kind 0..5 = Point, LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon).
pub fn golden_geometry(v: &Value) -> Geometry {
    let kind = v["kind"].as_i64().unwrap();
    let points: Vec<Coord> = v["points"].as_array().unwrap().iter().map(|s| seq(s)[0]).collect();
    let lines: Vec<Vec<Coord>> = v["lines"].as_array().unwrap().iter().map(seq).collect();
    let polys: Vec<Polygon> = v["polygons"]
        .as_array()
        .unwrap()
        .iter()
        .map(|p| Polygon::new(p.as_array().unwrap().iter().map(seq).collect()))
        .collect();
    match kind {
        0 => Geometry::Point(points[0]),
        1 => Geometry::LineString(lines.into_iter().next().unwrap()),
        2 => Geometry::Polygon(polys.into_iter().next().unwrap()),
        3 => Geometry::MultiPoint(points),
        4 => Geometry::MultiLineString(lines),
        5 => Geometry::MultiPolygon(polys),
        _ => panic!("bad kind"),
    }
}

pub fn golden_grid(v: &Value) -> GridSpec {
    let g = &v["grid"];
    GridSpec::new(
        g["xmin"].as_f64().unwrap(),
        g["ymin"].as_f64().unwrap(),
        g["xmax"].as_f64().unwrap(),
        g["ymax"].as_f64().unwrap(),
        g["ncol"].as_u64().unwrap() as u32,
        g["nrow"].as_u64().unwrap() as u32,
    )
}

pub fn golden_result(v: &Value) -> BurnResult {
    let i = |x: &Value| x.as_i64().unwrap() as i32;
    let f = |x: &Value| x.as_f64().unwrap() as f32;
    BurnResult {
        runs: v["runs"].as_array().unwrap().iter().map(|r| GridRun { row: i(&r[0]), col_start: i(&r[1]), col_end: i(&r[2]), id: i(&r[3]) }).collect(),
        edges: v["edges"].as_array().unwrap().iter().map(|r| GridEdge { row: i(&r[0]), col: i(&r[1]), fraction: f(&r[2]), id: i(&r[3]) }).collect(),
        lines: v["lines"].as_array().unwrap().iter().map(|r| GridLine { row: i(&r[0]), col: i(&r[1]), length: f(&r[2]), id: i(&r[3]) }).collect(),
        points: v["points"].as_array().unwrap().iter().map(|r| GridPoint { row: i(&r[0]), col: i(&r[1]), id: i(&r[2]) }).collect(),
        notes: Vec::new(),
    }
}

pub fn hex(s: &str) -> Vec<u8> {
    (0..s.len()).step_by(2).map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap()).collect()
}

pub fn square(x0: f64, y0: f64, x1: f64, y1: f64) -> Geometry {
    Geometry::Polygon(Polygon::new(vec![vec![
        Coord::new(x0, y0),
        Coord::new(x1, y0),
        Coord::new(x1, y1),
        Coord::new(x0, y1),
        Coord::new(x0, y0),
    ]]))
}

pub fn area(r: &BurnResult, g: &GridSpec) -> f64 {
    r.covered_cells() * g.dx() * g.dy()
}
