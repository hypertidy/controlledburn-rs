//! A hole's boundary cells inside the exterior's interior keep the
//! exterior's cover: coverage there is 1 minus the hole's fraction, and a
//! polygon with holes covers exactly its area.

use controlledburn::{burn, BurnOptions, Coord, Geometry, GridSpec, Polygon};

fn ring(v: &[(f64, f64)]) -> Vec<Coord> {
    v.iter().map(|&(x, y)| Coord::new(x, y)).collect()
}

fn area(v: &[(f64, f64)]) -> f64 {
    v.windows(2).map(|w| w[0].0 * w[1].1 - w[1].0 * w[0].1).sum::<f64>().abs() / 2.0
}

#[test]
fn holes_inside_the_exterior_cover_their_area() {
    let ext = [(0.5, 0.5), (9.5, 0.5), (9.5, 9.5), (0.5, 9.5), (0.5, 0.5)];
    let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
    let holes: [&[(f64, f64)]; 4] = [
        // across four cells, crossing row centres
        &[(4.2, 4.2), (5.8, 4.2), (5.8, 5.8), (4.2, 5.8), (4.2, 4.2)],
        // within one cell
        &[(4.2, 4.2), (4.8, 4.2), (4.8, 4.8), (4.2, 4.8), (4.2, 4.2)],
        // touching cells of a row without crossing its centre line
        &[(4.2, 4.6), (5.8, 4.6), (5.0, 5.4), (4.2, 4.6)],
        // the other winding
        &[(4.2, 4.2), (4.2, 5.8), (5.8, 5.8), (5.8, 4.2), (4.2, 4.2)],
    ];
    for hole in holes {
        let g = Geometry::Polygon(Polygon::new(vec![ring(&ext), ring(hole)]));
        let r = burn(&[g], &grid, BurnOptions::coverage()).unwrap();
        let want = area(&ext) - area(hole);
        assert!((r.covered_cells() - want).abs() < 1e-5, "covered {} of {want}", r.covered_cells());
        // the cell at (row 5, col 4) holds part of every hole: partly covered
        let e = r.edges.iter().find(|e| e.row == 5 && e.col == 4).expect("an edge cell at (5, 4)");
        assert!(e.fraction > 0.0 && e.fraction < 1.0);
    }
}

#[test]
fn two_holes_and_an_island_of_cells() {
    let ext = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)];
    let h1 = [(1.3, 1.3), (3.7, 1.3), (3.7, 3.7), (1.3, 3.7), (1.3, 1.3)];
    let h2 = [(6.1, 6.6), (8.9, 6.6), (7.5, 8.8), (6.1, 6.6)];
    let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
    let g = Geometry::Polygon(Polygon::new(vec![ring(&ext), ring(&h1), ring(&h2)]));
    let r = burn(&[g], &grid, BurnOptions::coverage()).unwrap();
    let want = area(&ext) - area(&h1) - area(&h2);
    assert!((r.covered_cells() - want).abs() < 1e-5, "covered {} of {want}", r.covered_cells());
}

#[test]
fn many_holes_cover_their_area() {
    // 20 x 20 diamond holes, none aligned with the grid. Also guards the
    // per-row cost: the enclosing-ring sum must not loop over every ring
    // for every boundary cell.
    let ext = [(0.0, 0.0), (100.0, 0.0), (100.0, 100.0), (0.0, 100.0), (0.0, 0.0)];
    let mut rings = vec![ring(&ext)];
    let mut want = area(&ext);
    for i in 0..20 {
        for j in 0..20 {
            let (cx, cy, r) = (2.5 + 5.0 * i as f64 + 0.013, 2.5 + 5.0 * j as f64 + 0.029, 1.7);
            let h = [(cx - r, cy), (cx, cy - r), (cx + r, cy), (cx, cy + r), (cx - r, cy)];
            want -= area(&h);
            rings.push(ring(&h));
        }
    }
    let grid = GridSpec::new(0.0, 0.0, 100.0, 100.0, 333, 333);
    let r = burn(&[Geometry::Polygon(Polygon::new(rings))], &grid, BurnOptions::coverage()).unwrap();
    let got = r.covered_cells() * grid.dx() * grid.dy();
    assert!((got - want).abs() < 1e-3, "covered {got} of {want}");
}
