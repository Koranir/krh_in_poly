#![cfg_attr(feature = "no_std", no_std)]
#![warn(clippy::all, clippy::pedantic, clippy::nursery)]
#![allow(clippy::float_cmp)]

#[cfg(feature = "no_std")]
#[allow(clippy::inline_always)]
mod math {
    #[inline(always)]
    pub fn sqrt(a: f32) -> f32 {
        libm::sqrtf(a)
    }

    #[inline(always)]
    pub fn mul_add(a: f32, mul: f32, add: f32) -> f32 {
        a * mul + add
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
    pub fn mul_add(a: f32, mul: f32, add: f32) -> f32 {
        a.mul_add(mul, add)
    }
}

use math::{mul_add, sqrt};

// #[must_use]
// #[allow(clippy::similar_names)]
// fn cubic_bezier_wind(
//     (px, py): (f32, f32),
//     (qax, qay): (f32, f32),
//     (qbx, qby): (f32, f32),
//     (qcx, qcy): (f32, f32),
//     (qdx, qdy): (f32, f32),
// ) -> i32 {
//     // solve p = ((1-t)^3)a + (3t(1-t)^2)b + (3t^2)(1-t)c + (t^3)d for t

//     todo!()
// }

#[must_use]
#[allow(clippy::similar_names)]
pub fn quadratic_bezier_wind(
    (px, py): (f32, f32),
    (qax, qay): (f32, f32),
    (qbx, qby): (f32, f32),
    (qcx, qcy): (f32, f32),
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

#[must_use]
/// Check if a point is in a polygon.
pub fn in_polygon((px, py): (f32, f32), poly: &[(f32, f32)]) -> bool {
    // Iterator of the points in the polygon.
    let poly_iter = poly.iter().copied().cycle();

    let winding_number: i32 = poly_iter
        .clone()
        // The starts of the edges.
        .take(poly.len())
        // The ends of the edges.
        .zip(poly_iter.skip(1).take(poly.len()))
        // Only calculate a winding number for edges that are in the y and x bounds.
        .map(|((x0, y0), (x1, y1))| {
            // A point `py` is in the y bounds if it is in `y0..=y1`.
            let in_y = (py <= y0) == (py >= y1);
            // Make sure divide by zeros don't break things.
            let my = if y0 == y1 { 1.0 } else { (y0 - py) / (y0 - y1) };
            // Linear segment, can linearly interpolate between x0 and x1 with the percentage between y0 and y1 py is to get the x value of the edge at py.
            let vx = x1 - x0;
            let in_x = px < mul_add(vx, my, x0); // (x0 + vx * my);

            // Don't change the winding number if we're not in bounds of a line.
            if !(in_y && in_x) {
                return 0;
            }
            // Special case if we're directly on an edge/point.
            if in_y && in_x && py == y0 || px == x0 {
                return 1;
            }
            // Winding number, when we approach a line from the left ( -> ):
            // If the line's last point is top-right or bottom-left ( / ), increasing
            // If the line's last point is top-left or bottom-right ( \ ), decreasing
            match (x0 > x1, y0 > y1) {
                (true, true) | (false, false) => 1,
                (true, false) | (false, true) => -1,
            }
        })
        .sum();

    // If the winding number is zero, we are outside the polygon.
    winding_number != 0
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
}
