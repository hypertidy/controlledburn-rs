//! A cell whose only edges lie along its walls, one crossing the row centre
//! up and the other down, has zero coverage and zero net winding. It is
//! still a boundary cell and must not be filled as part of an interior run.

use controlledburn::{burn, BurnOptions, Coord, Geometry, GridSpec, Polygon};

#[test]
fn a_cell_with_edges_only_along_its_walls_stays_a_boundary_cell() {
    // A lattice-aligned notch, from a raster-derived waterbody (DEA
    // Waterbodies uid q9g5b3eyr_v1). Cell [4, 5] x [2, 3] has edges along
    // both walls; it used to be filled, 3.5 covered cells instead of 2.5.
    let v = [(4.0, 2.0), (2.0, 7.0), (6.0, 2.0), (5.0, 2.0), (5.0, 3.0), (3.0, 5.5), (4.0, 2.5), (4.0, 2.0)];
    let ring: Vec<Coord> = v.iter().map(|&(x, y)| Coord::new(x, y)).collect();
    let g = Geometry::Polygon(Polygon::new(vec![ring]));
    let grid = GridSpec::new(0.0, 0.0, 10.0, 10.0, 10, 10);
    let r = burn(&[g], &grid, BurnOptions::coverage()).unwrap();
    assert!((r.covered_cells() - 2.5).abs() < 1e-5, "covered {} of 2.5", r.covered_cells());
    // the notch cell is row 7, col 4 (0-based): outside
    assert!(!r.runs.iter().any(|run| run.row == 7 && run.col_start <= 4 && run.col_end > 4));
    assert!(!r.edges.iter().any(|e| e.row == 7 && e.col == 4));
}

#[test]
fn a_zero_width_spike_into_a_covered_cell_keeps_it_covered() {
    // Found by the property suite: the ring runs in through the top-left
    // corner of cell [4.8, 5.1] x [1.2, 1.5] to its centre and straight back
    // out the same way. That traversal leaves where it entered and encloses
    // nothing; the cell is inside the polygon and must stay covered.
    let v = [
        (2.2955346720684804, -0.28831957426678434),
        (4.2871984006926445, -0.6968703145658899),
        (6.183801073287196, -1.4711665488619694),
        (5.24135385306532, 1.4474746273177665),
        (5.275730794077829, 1.5426958381891471),
        (5.824891418164306, 2.019485421020327),
        (5.151353853065319, 1.696766993828962),
        (3.5999999999999996, 2.6999999999999997),
        (4.95, 1.3499999999999996),
        (4.8, 1.4999999999999996),
        (2.2955346720684804, -0.28831957426678434),
    ];
    let ring: Vec<Coord> = v.iter().map(|&(x, y)| Coord::new(x, y)).collect();
    let g = Geometry::Polygon(Polygon::new(vec![ring]));
    let grid = GridSpec::new(0.0, 0.0, 6.0, 5.3999999999999995, 20, 18);
    let r = burn(&[g], &grid, BurnOptions::coverage()).unwrap();
    // the spike cell is row 13, col 16 (0-based): fully covered
    assert!(r.runs.iter().any(|run| run.row == 13 && run.col_start <= 16 && run.col_end > 16));
    // area of the polygon clipped to the grid (cells of 0.3 x 0.3)
    assert!((r.covered_cells() * 0.09 - 3.026275470389206).abs() < 1e-5, "covered {}", r.covered_cells() * 0.09);
}
