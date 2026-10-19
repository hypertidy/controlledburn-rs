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
