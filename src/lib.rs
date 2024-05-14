#![no_std]
#![warn(clippy::all, clippy::pedantic, clippy::nursery)]
#![allow(clippy::float_cmp)]

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
        .filter_map(|((x0, y0), (x1, y1))| {
            // A point `py` is in the y bounds if it is in `y0..=y1`.
            let in_y = (py <= y0) == (py >= y1);
            // Make sure divide by zeros don't break things.
            let my = if y0 == y1 { 1.0 } else { (y0 - py) / (y0 - y1) };
            // Linear segment, can linearly interpolate between x0 and x1 with the percentage between y0 and y1 py is to get the x value of the edge at py.
            let vx = x1 - x0;
            let in_x = px < (x0 + vx * my);
            // Special case if we're directly on an edge/point.
            if in_y && in_x && py == y0 || px == x0 {
                return Some(1);
            }
            // Winding number, when we approach a line from the left ( -> ):
            // If the line's last point is top-right or bottom-left ( / ), increasing
            // If the line's last point is top-left or bottom-right ( \ ), decreasing
            (in_y && in_x).then_some(match (x0 > x1, y0 > y1) {
                (true, true) | (false, false) => 1,
                (true, false) | (false, true) => -1,
            })
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
}
