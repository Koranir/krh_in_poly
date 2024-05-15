#![cfg_attr(feature = "no_std", no_std)]
#![warn(clippy::all, clippy::pedantic, clippy::nursery)]
#![allow(clippy::float_cmp)]

use math::{arccos, cbrt, cos, mul_add, sqrt};

#[cfg(feature = "no_std")]
#[allow(clippy::inline_always)]
mod math {
    #[inline(always)]
    pub fn sqrt(a: f32) -> f32 {
        libm::sqrtf(a)
    }

    #[inline(always)]
    pub fn cbrt(a: f32) -> f32 {
        libm::cbrtf(a)
    }

    #[inline(always)]
    pub fn mul_add(a: f32, mul: f32, add: f32) -> f32 {
        a * mul + add
    }

    #[inline(always)]
    pub fn arccos(a: f32) -> f32 {
        libm::acosf(a)
    }

    #[inline(always)]
    pub fn cos(a: f32) -> f32 {
        libm::cosf(a)
    }
}

#[cfg(feature = "std")]
#[allow(clippy::inline_always)]
mod math {
    #[inline(always)]
    pub fn sqrt(a: f32) -> f32 {
        a.sqrt()
    }

    #[inline(always)]
    pub fn cbrt(a: f32) -> f32 {
        a.cbrt()
    }

    #[inline(always)]
    pub fn mul_add(a: f32, mul: f32, add: f32) -> f32 {
        a.mul_add(mul, add)
    }

    #[inline(always)]
    pub fn arccos(a: f32) -> f32 {
        a.acos()
    }

    #[inline(always)]
    pub fn cos(a: f32) -> f32 {
        a.cos()
    }
}

/// A point on the cartesian plane.
/// (x, y)
pub type Point = (f32, f32);
/// A line segment.
/// (start, end)
pub type Line = (Point, Point);
/// A quadratic bezier.
/// (start, control, end)
pub type QuadBezier = (Point, Point, Point);
/// A cubic bezier.
/// (start, start control, end control, end)
pub type CubeBezier = (Point, Point, Point, Point);

#[must_use]
#[allow(clippy::similar_names, clippy::many_single_char_names)]
#[inline]
pub fn cubic_bezier_wind(
    (px, py): Point,
    (ax, ay): Point,
    (bx, by): Point,
    (cx, cy): Point,
    (dx, dy): Point,
) -> i32 {
    let cubic_bezier_x_at = |t| {
        let m1t = 1.0 - t;
        let m1t2 = m1t * m1t;
        let t2 = t * t;
        let a = ax * m1t2 * m1t;
        let b = bx * 3.0 * t * m1t2;
        let c = 3.0 * t2 * m1t * cx;
        let d = t2 * t * dx;

        a + b + c + d
    };

    let cubic_bezier_ddx_at = |t| {
        let m1t = 1.0 - t;
        let m1t2 = m1t * m1t;
        let t2 = t * t;
        let tm1t2 = 3.0 * m1t2;
        let sm1tt = 6.0 * m1t * t;
        let tt2 = 3.0 * t2;
        let ay = tm1t2 * (by - ay);
        let ax = tm1t2 * (bx - ax);
        let by = sm1tt * (cy - by);
        let bx = sm1tt * (cx - bx);
        let cy = tt2 * (dy - cy);
        let cx = tt2 * (dx - cx);

        let dy = ay + by + cy;
        let dx = ax + bx + cx;

        dy / dx
    };

    // let a = -ay + (3.0 * (by - cy)) + dy;
    let a = mul_add(by - cy, 3.0, dy) - ay;
    // let b = 3.0 * (ay - 2.0 * by + cy);
    let b = 3.0 * by.mul_add(-2.0, ay + cy);
    let c = 3.0 * (by - ay);
    let d = ay - py;

    let a2 = a * a;
    let b2 = b * b;
    // let p = ((3.0 * a * c) - b2) / (3.0 * a2);
    let p = mul_add(3.0 * a, c, -b2) / (3.0 * a2);
    let qa227 = 27.0 * a2;
    // let qn = (2.0 * b * b2) - (9.0 * a * b * c) + (qa227 * d);
    let qn = mul_add(b2, b * 2.0, mul_add(qa227, d, -9.0 * a * b * c));
    let q = qn / (qa227 * a);

    let correction = b / (-3.0 * a);

    let p2 = p * p;
    let p3 = p2 * p;
    let q2 = q * q;

    // let discrim = (-4.0 * p3) + (-27.0 * q2);
    let discrim = mul_add(p3, -4.0, q2 * -27.0);

    let (root, others) = if discrim > 0.0 {
        let tk0 = 2.0 * sqrt(-p / 3.0);
        let tk1 = (3.0 * q) / (2.0 * p);
        let tk1 = tk1 * sqrt(-3.0 / p);
        let tk1 = (1.0 / 3.0) * arccos(tk1);

        let root1 = tk0 * cos(tk1);
        // let root2 = tk0 * cos(tk1 - core::f32::consts::TAU * 2.0 * (1.0 / 3.0));
        #[allow(clippy::unreadable_literal)]
        let root2 = tk0 * cos(tk1 - 4.1887902);
        // let root3 = tk0 * cos(tk1 - core::f32::consts::TAU * 2.0 * (2.0 / 3.0));
        let root3 = tk0 * cos(tk1 - 8.37758);

        (
            root1 + correction,
            Some((root2 + correction, root3 + correction)),
        )
    } else {
        let tc0 = sqrt(q2 / 4.0 + p3 / 27.0);
        let mq2 = q / -2.0;

        let tc1 = cbrt(mq2 + tc0);
        let tc2 = cbrt(mq2 - tc0);

        (tc1 + tc2 + correction, None)
    };

    let mut wind = 0;
    if (0.0..=1.0).contains(&root) {
        let qtx0 = cubic_bezier_x_at(root);
        if px <= qtx0 {
            let qddx0 = cubic_bezier_ddx_at(root);
            wind += if qddx0 < 0.0 { 1 } else { -1 };
        }
    }
    if let Some((r1, r2)) = others {
        if (0.0..=1.0).contains(&r1) {
            let qtx0 = cubic_bezier_x_at(r1);
            if px <= qtx0 {
                let qddx0 = cubic_bezier_ddx_at(r1);
                wind += if qddx0 < 0.0 { 1 } else { -1 };
            }
        }
        if (0.0..=1.0).contains(&r2) {
            let qtx0 = cubic_bezier_x_at(r2);
            if px <= qtx0 {
                let qddx0 = cubic_bezier_ddx_at(r2);
                wind += if qddx0 < 0.0 { 1 } else { -1 };
            }
        }
    }

    wind
}

#[must_use]
#[inline]
#[allow(clippy::similar_names)]
pub fn quadratic_bezier_wind(
    (px, py): Point,
    (qax, qay): Point,
    (qbx, qby): Point,
    (qcx, qcy): Point,
) -> i32 {
    let det = mul_add(
        qay,
        py,
        mul_add(qby, qby, mul_add(qay, qcy, -2.0 * qby * py)),
    ); // (qay * py) + (qby * qby) + (qay * qcy) - (2.0 * qby * py);
    if det < 0.0 {
        return 0;
    }

    let ab = qay - qby;
    let a2bc = mul_add(qby, -2.0, qay) + qcy; //  qay - (2.0 * qby) + qcy;

    let t0 = (sqrt(det) + ab) / a2bc;
    let t1 = (-sqrt(det) + ab) / a2bc;

    let quad_x = |t| {
        let tm1 = 1.0 - t;
        // (tm1 * tm1 * qax) + (2.0 * t * tm1 * qbx) + (t * t * qcx)
        mul_add(t * t, qcx, mul_add(tm1 * tm1, qax, 2.0 * t * tm1 * qbx))
    };
    let quad_ddx = |t| {
        let trcp1 = 1.0 / t - 1.0;

        // let dy = trcp1 * (qby - qay) + (qcy - qby);
        // let dx = trcp1 * (qbx - qax) + (qcx - qbx);
        let dy = mul_add(trcp1, qby - qay, qcy - qby);
        let dx = mul_add(trcp1, qbx - qax, qcx - qbx);

        dy / dx
    };

    let mut wind = 0;
    if (0.0..=1.0).contains(&t0) {
        let qtx0 = quad_x(t0);
        if px <= qtx0 {
            let qddx0 = quad_ddx(t0);
            wind += if qddx0 < 0.0 { 1 } else { -1 };
        }
    }
    if (0.0..=1.0).contains(&t1) {
        let qtx1 = quad_x(t1);

        if px <= qtx1 {
            let qddx1 = quad_ddx(qtx1);
            wind += if qddx1 < 0.0 { 1 } else { -1 };
        }
    }

    wind
}

#[inline]
#[must_use]
pub fn line_segment_wind((px, py): Point, (x0, y0): Point, (x1, y1): Point) -> i32 {
    // A point `py` is in the y bounds if it is in `y0..=y1`.
    let in_y = (py <= y0) == (py >= y1);
    // Make sure divide by zeros don't break things.
    let my = if y0 == y1 { 1.0 } else { (y0 - py) / (y0 - y1) };
    // Linear segment, can linearly interpolate between x0 and x1 with the percentage between y0 and y1 py is to get the x value of the edge at py.
    let vx = x1 - x0;
    let in_x = px < mul_add(vx, my, x0); // (x0 + vx * my);

    // Special case if we're directly on an edge/point.
    if in_y && in_x && py == y0 || px == x0 {
        return 1;
    }
    // Don't change the winding number if we're not in bounds of a line.
    if !(in_y && in_x) {
        return 0;
    }
    // Winding number, when we approach a line from the left ( -> ):
    // If the line's last point is top-right or bottom-left ( / ), increasing
    // If the line's last point is top-left or bottom-right ( \ ), decreasing
    match (x0 > x1, y0 > y1) {
        (true, true) | (false, false) => 1,
        (true, false) | (false, true) => -1,
    }
}

#[must_use]
/// Check if a point is in a polygon.
pub fn in_polygon((px, py): Point, poly: &[Point]) -> bool {
    // Iterator of the points in the polygon.
    let poly_iter = poly.iter().copied().cycle();

    let winding_number: i32 = poly_iter
        .clone()
        // The starts of the edges.
        .take(poly.len())
        // The ends of the edges.
        .zip(poly_iter.skip(1).take(poly.len()))
        // Only calculate a winding number for edges that are in the y and x bounds.
        .map(|(start, end)| line_segment_wind((px, py), start, end))
        .sum();

    // If the winding number is zero, we are outside the polygon.
    winding_number != 0
}

#[must_use]
pub fn in_poly_with_curves(
    (px, py): (f32, f32),
    lines: &[Line],
    quadratics: &[QuadBezier],
    cubics: &[CubeBezier],
) -> bool {
    let line_wind: i32 = lines
        .iter()
        .map(|f| line_segment_wind((px, py), f.0, f.1))
        .sum();
    let quad_wind: i32 = quadratics
        .iter()
        .map(|f| quadratic_bezier_wind((px, py), f.0, f.1, f.2))
        .sum();
    let cube_wind: i32 = cubics
        .iter()
        .map(|f| cubic_bezier_wind((px, py), f.0, f.1, f.2, f.3))
        .sum();

    let wind = line_wind + quad_wind + cube_wind;

    wind != 0
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_point_inside_polygon() {
        assert!(in_polygon(
            (0.5, 0.5),
            &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        ));
    }

    #[test]
    fn test_point_inside_polygon_rev() {
        assert!(in_polygon(
            (0.5, 0.5),
            &[(0.0, 1.0), (1.0, 1.0), (1.0, 0.0), (0.0, 0.0)]
        ));
    }

    #[test]
    fn test_point_outside_polygon() {
        assert!(!in_polygon(
            (2.0, 2.0),
            &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        ));
    }

    #[test]
    fn test_point_on_edge() {
        assert!(in_polygon(
            (0.0, 0.0),
            &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        ));
    }

    #[test]
    fn test_point_inside_complex_polygon() {
        assert!(in_polygon(
            (1.0, 1.0),
            &[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (1.0, 2.0), (0.0, 1.0)]
        ));
    }

    #[test]
    fn test_point_outside_complex_polygon() {
        assert!(!in_polygon(
            (2.5, 2.5),
            &[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (1.0, 2.0), (0.0, 1.0)]
        ));
    }

    #[test]
    fn test_point_on_vertex() {
        assert!(in_polygon(
            (1.0, 0.0),
            &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        ));
    }

    #[test]
    fn test_point_outside_extension_of_edge() {
        assert!(!in_polygon(
            (-1.0, 0.0),
            &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        ));
    }

    #[test]
    fn test_point_inside_concave_polygon() {
        assert!(in_polygon(
            (0.5, 0.5),
            &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.5, 0.5), (0.0, 1.0)]
        ));
    }
    #[test]
    fn test_point_outside_concave_polygon() {
        assert!(!in_polygon(
            (0.25, -0.25),
            &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.5, 0.5), (0.0, 1.0)]
        ));
    }

    #[test]
    fn test_point_on_boundary_of_concave_polygon() {
        assert!(in_polygon(
            (0.0, 0.0),
            &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.5, 0.5), (0.0, 1.0)]
        ));
    }

    #[test]
    fn test_point_outside_large_polygon() {
        assert!(!in_polygon(
            (-5.0, 5.0),
            &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
        ));
    }

    #[test]
    fn test_point_inside_large_polygon() {
        assert!(in_polygon(
            (5.0, 5.0),
            &[
                (0.0, 0.0),
                (10.0, 0.0),
                (10.0, 10.0),
                (5.0, 5.0),
                (0.0, 2.0)
            ]
        ));
    }

    #[test]
    fn test_point_on_edge_large_polygon() {
        assert!(in_polygon(
            (0.0, 0.0),
            &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
        ));
    }

    #[test]
    fn test_point_outside_polygon_with_negative_coordinates() {
        assert!(!in_polygon(
            (-1.0, -1.0),
            &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        ));
    }

    #[test]
    fn test_point_inside_polygon_with_negative_coordinates() {
        assert!(in_polygon(
            (-0.5, -0.5),
            &[(0.0, 0.0), (-1.0, 0.0), (-1.0, -1.0), (0.0, -1.0)]
        ));
    }
    #[test]
    fn test_point_on_boundary() {
        assert!(in_polygon(
            (1.0, 0.0),
            &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        )); // On horizontal edge
        assert!(in_polygon(
            (1.0, 1.0),
            &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        )); // On vertical edge
    }

    #[test]
    fn test_point_outside_polygon_close_to_edge() {
        assert!(!in_polygon(
            (-0.1, -0.1),
            &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        )); // Close to corner
        assert!(!in_polygon(
            (0.5, 1.1),
            &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        )); // Close to edge
    }

    #[test]
    fn test_point_on_boundary_of_complex_polygon() {
        assert!(in_polygon(
            (0.5, 0.0),
            &[(0.0, 0.0), (1.0, 0.0), (2.0, 1.0), (1.0, 2.0), (0.0, 2.0)]
        )); // On horizontal edge
        assert!(in_polygon(
            (1.0, 1.5),
            &[(0.0, 0.0), (1.0, 0.0), (2.0, 1.0), (1.0, 2.0), (0.0, 2.0)]
        )); // On vertical edge
    }

    #[test]
    fn test_point_outside_complex_polygon_close_to_edge() {
        assert!(!in_polygon(
            (-0.1, -0.1),
            &[(0.0, 0.0), (1.0, 0.0), (2.0, 1.0), (1.0, 2.0), (0.0, 2.0)]
        )); // Close to corner
        assert!(!in_polygon(
            (0.5, 2.1),
            &[(0.0, 0.0), (1.0, 0.0), (2.0, 1.0), (1.0, 2.0), (0.0, 2.0)]
        )); // Close to edge
    }

    #[test]
    fn test_point_on_vertex_of_complex_polygon() {
        assert!(in_polygon(
            (0.0, 0.0),
            &[(0.0, 0.0), (1.0, 0.0), (2.0, 1.0), (1.0, 2.0), (0.0, 2.0)]
        )); // On vertex
        assert!(in_polygon(
            (2.0, 1.0),
            &[(0.0, 0.0), (1.0, 0.0), (2.0, 1.0), (1.0, 2.0), (0.0, 2.0)]
        )); // On vertex
    }

    #[test]
    fn test_point_inside_polygon_close_to_edge() {
        assert!(in_polygon(
            (0.01, 0.01),
            &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        )); // Close to corner
        assert!(in_polygon(
            (0.5, 0.99),
            &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        )); // Close to edge
    }

    #[test]
    fn test_point_outside_polygon_far_from_edge() {
        assert!(!in_polygon(
            (-1.0, -1.0),
            &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        )); // Far from any edge
        assert!(!in_polygon(
            (10.0, 10.0),
            &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        )); // Far from any edge
    }
    #[test]
    fn test_empty_polygon() {
        // An empty polygon should always return false
        assert!(!in_polygon((0.0, 0.0), &[]));
    }

    #[test]
    fn test_single_point_polygon() {
        // A polygon with a single point should always return false for any point not equal to that point
        assert!(!in_polygon((0.0, 0.0), &[(1.0, 1.0)]));
        assert!(!in_polygon((2.0, 2.0), &[(1.0, 1.0)]));
        assert!(!in_polygon((-1.0, -1.0), &[(0.0, 0.0)]));
    }

    #[test]
    fn test_horizontal_line_polygon() {
        // A horizontal line polygon should return true for points lying on the line and false otherwise
        assert!(in_polygon((0.5, 0.0), &[(0.0, 0.0), (1.0, 0.0)])); // On the line
        assert!(!in_polygon((0.5, 1.0), &[(0.0, 0.0), (1.0, 0.0)])); // Above the line
        assert!(!in_polygon((0.5, -1.0), &[(0.0, 0.0), (1.0, 0.0)])); // Below the line
    }

    #[test]
    fn test_vertical_line_polygon() {
        // A vertical line polygon should return true for points lying on the line and false otherwise
        assert!(in_polygon((0.0, 0.5), &[(0.0, 0.0), (0.0, 1.0)])); // On the line
        assert!(!in_polygon((1.0, 0.5), &[(0.0, 0.0), (0.0, 1.0)])); // Right to the line
        assert!(!in_polygon((-1.0, 0.5), &[(0.0, 0.0), (0.0, 1.0)])); // Left to the line
    }

    #[test]
    fn test_large_polygon() {
        // Test with a large polygon with many vertices
        let polygon = [
            (0.0, 0.0),
            (1.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (1.0, 1.0),
            (1.0, 2.0),
            (0.0, 2.0),
        ];
        assert!(in_polygon((0.5, 1.5), &polygon)); // Inside the polygon
        assert!(!in_polygon((3.0, 1.0), &polygon)); // Outside the polygon
    }

    #[test]
    fn test_negative_coordinates_polygon() {
        // Test with a polygon containing negative coordinates
        let polygon = [(-1.0, -1.0), (-1.0, 1.0), (1.0, 1.0), (1.0, -1.0)];
        assert!(in_polygon((0.0, 0.0), &polygon)); // Inside the polygon
        assert!(!in_polygon((2.0, 2.0), &polygon)); // Outside the polygon
    }

    #[test]
    fn test_bezier_curve_wind() {
        assert_eq!(
            quadratic_bezier_wind((-3.0, 2.0), (-1.0, 4.0), (-2.0, -3.0), (4.0, 1.0)),
            -1
        );
    }

    #[test]
    fn test_cubic_curve_wind() {
        assert_eq!(
            cubic_bezier_wind(
                (-20.0, 0.0),
                (31.4, -25.0),
                (24.6, 9.9),
                (-12.5, -4.6),
                (-20.0, -29.0)
            ),
            0
        );
        assert_eq!(
            cubic_bezier_wind(
                (0.0, -20.0),
                (31.4, -25.0),
                (24.6, 9.9),
                (-18.0, 24.6),
                (-20.0, -29.0)
            ),
            1
        );
    }
}
