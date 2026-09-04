//! Point rasterization: one record per point inside the grid.

use crate::geometry::Coord;
use crate::grid::Grid;
use crate::output::GridPoint;

/// A point is either in a cell or it is not, so points carry no weight.
/// All four extent edges are inclusive: a point exactly on `xmax` or
/// `ymin` goes to the last column/row. Points outside are dropped.
pub(crate) fn process_point(pt: &Coord, full_grid: &Grid<0>, point_id: i32, out: &mut Vec<GridPoint>) {
    if !pt.is_finite() {
        return; // POINT EMPTY convention
    }
    let ext = full_grid.extent();
    if pt.x < ext.xmin || pt.x > ext.xmax || pt.y < ext.ymin || pt.y > ext.ymax {
        return;
    }
    let r = full_grid.get_row(pt.y);
    let c = full_grid.get_column(pt.x);
    out.push(GridPoint { row: r as i32, col: c as i32, id: point_id });
}
