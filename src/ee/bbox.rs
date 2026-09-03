// Derived from exactextract box.h / box.cpp / crossing.h
// Copyright (c) 2018-2019 ISciences, LLC. Apache License 2.0.

use super::side::Side;
use crate::geometry::Coord;

/// An axis-aligned box.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BBox {
    pub xmin: f64,
    pub ymin: f64,
    pub xmax: f64,
    pub ymax: f64,
}

/// Where a segment leaves a box: the side and the exit coordinate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Crossing {
    pub side: Side,
    pub coord: Coord,
}

#[inline]
fn clamp(x: f64, low: f64, high: f64) -> f64 {
    // std::min(std::max(x, low), high) with C++ semantics
    let m = if low > x { low } else { x };
    if high < m {
        high
    } else {
        m
    }
}

impl BBox {
    #[inline]
    pub const fn new(xmin: f64, ymin: f64, xmax: f64, ymax: f64) -> Self {
        BBox { xmin, ymin, xmax, ymax }
    }

    #[inline]
    pub const fn make_empty() -> Self {
        BBox::new(0.0, 0.0, 0.0, 0.0)
    }

    #[inline]
    pub fn width(&self) -> f64 {
        self.xmax - self.xmin
    }

    #[inline]
    pub fn height(&self) -> f64 {
        self.ymax - self.ymin
    }

    #[inline]
    pub fn area(&self) -> f64 {
        self.width() * self.height()
    }

    #[inline]
    pub fn perimeter(&self) -> f64 {
        2.0 * self.width() + 2.0 * self.height()
    }

    pub fn intersects(&self, other: &BBox) -> bool {
        !(other.ymin > self.ymax || other.ymax < self.ymin || other.xmin > self.xmax || other.xmax < self.xmin)
    }

    pub fn intersection(&self, other: &BBox) -> BBox {
        BBox::new(
            self.xmin.max(other.xmin),
            self.ymin.max(other.ymin),
            self.xmax.min(other.xmax),
            self.ymax.min(other.ymax),
        )
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.xmin >= self.xmax || self.ymin >= self.ymax
    }

    #[inline]
    pub fn contains(&self, c: &Coord) -> bool {
        c.x >= self.xmin && c.x <= self.xmax && c.y >= self.ymin && c.y <= self.ymax
    }

    #[inline]
    pub fn strictly_contains(&self, c: &Coord) -> bool {
        c.x > self.xmin && c.x < self.xmax && c.y > self.ymin && c.y < self.ymax
    }

    /// Which side a boundary coordinate lies on, with the exactextract
    /// priority Left, Right, Bottom, Top (a corner reports Left or Right).
    /// `None` if the coordinate is not on the boundary.
    pub fn side(&self, c: &Coord) -> Option<Side> {
        if c.x == self.xmin {
            Some(Side::Left)
        } else if c.x == self.xmax {
            Some(Side::Right)
        } else if c.y == self.ymin {
            Some(Side::Bottom)
        } else if c.y == self.ymax {
            Some(Side::Top)
        } else {
            None
        }
    }

    /// Exit crossing of the segment `c1 -> c2`, where `c2` lies outside
    /// the box. Returns `None` in the cases the C++ throws "Never get
    /// here" (an axis-parallel segment whose target is inside the box).
    pub fn crossing(&self, c1: &Coord, c2: &Coord) -> Option<Crossing> {
        // vertical line
        if c1.x == c2.x {
            return if c2.y >= self.ymax {
                Some(Crossing { side: Side::Top, coord: Coord::new(c1.x, self.ymax) })
            } else if c2.y <= self.ymin {
                Some(Crossing { side: Side::Bottom, coord: Coord::new(c1.x, self.ymin) })
            } else {
                None
            };
        }
        // horizontal line
        if c1.y == c2.y {
            return if c2.x >= self.xmax {
                Some(Crossing { side: Side::Right, coord: Coord::new(self.xmax, c1.y) })
            } else if c2.x <= self.xmin {
                Some(Crossing { side: Side::Left, coord: Coord::new(self.xmin, c1.y) })
            } else {
                None
            };
        }

        let m = ((c2.y - c1.y) / (c2.x - c1.x)).abs();
        let up = c2.y > c1.y;
        let right = c2.x > c1.x;

        let cr = if up {
            if right {
                // 1st quadrant
                let y2 = c1.y + m * (self.xmax - c1.x);
                if y2 < self.ymax {
                    Crossing { side: Side::Right, coord: Coord::new(self.xmax, clamp(y2, self.ymin, self.ymax)) }
                } else {
                    let x2 = c1.x + (self.ymax - c1.y) / m;
                    Crossing { side: Side::Top, coord: Coord::new(clamp(x2, self.xmin, self.xmax), self.ymax) }
                }
            } else {
                // 2nd quadrant
                let y2 = c1.y + m * (c1.x - self.xmin);
                if y2 < self.ymax {
                    Crossing { side: Side::Left, coord: Coord::new(self.xmin, clamp(y2, self.ymin, self.ymax)) }
                } else {
                    let x2 = c1.x - (self.ymax - c1.y) / m;
                    Crossing { side: Side::Top, coord: Coord::new(clamp(x2, self.xmin, self.xmax), self.ymax) }
                }
            }
        } else if right {
            // 4th quadrant
            let y2 = c1.y - m * (self.xmax - c1.x);
            if y2 > self.ymin {
                Crossing { side: Side::Right, coord: Coord::new(self.xmax, clamp(y2, self.ymin, self.ymax)) }
            } else {
                let x2 = c1.x + (c1.y - self.ymin) / m;
                Crossing { side: Side::Bottom, coord: Coord::new(clamp(x2, self.xmin, self.xmax), self.ymin) }
            }
        } else {
            // 3rd quadrant
            let y2 = c1.y - m * (c1.x - self.xmin);
            if y2 > self.ymin {
                Crossing { side: Side::Left, coord: Coord::new(self.xmin, clamp(y2, self.ymin, self.ymax)) }
            } else {
                let x2 = c1.x - (c1.y - self.ymin) / m;
                Crossing { side: Side::Bottom, coord: Coord::new(clamp(x2, self.xmin, self.xmax), self.ymin) }
            }
        };
        Some(cr)
    }

    #[inline]
    pub fn corners_ccw_from_bottom_left(&self) -> [Coord; 4] {
        [
            Coord::new(self.xmin, self.ymin),
            Coord::new(self.xmin, self.ymax),
            Coord::new(self.xmax, self.ymax),
            Coord::new(self.xmax, self.ymin),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crossing_quadrants() {
        let b = BBox::new(0., 0., 1., 1.);
        let c = b.crossing(&Coord::new(0.5, 0.5), &Coord::new(2.0, 1.0)).unwrap();
        assert_eq!(c.side, Side::Right);
        assert_eq!(c.coord, Coord::new(1.0, 0.5 + (1.0 / 3.0) * 0.5));
        let c = b.crossing(&Coord::new(0.5, 0.5), &Coord::new(0.5, 3.0)).unwrap();
        assert_eq!((c.side, c.coord), (Side::Top, Coord::new(0.5, 1.0)));
        let c = b.crossing(&Coord::new(0.5, 0.5), &Coord::new(-1.0, -1.0)).unwrap();
        assert_eq!((c.side, c.coord), (Side::Bottom, Coord::new(0.0, 0.0)));
        assert!(b.crossing(&Coord::new(0.5, 0.5), &Coord::new(0.5, 0.7)).is_none());
    }

    #[test]
    fn side_priority() {
        let b = BBox::new(0., 0., 1., 1.);
        assert_eq!(b.side(&Coord::new(0., 1.)), Some(Side::Left));
        assert_eq!(b.side(&Coord::new(1., 0.)), Some(Side::Right));
        assert_eq!(b.side(&Coord::new(0.5, 0.)), Some(Side::Bottom));
        assert_eq!(b.side(&Coord::new(0.5, 0.5)), None);
    }
}
