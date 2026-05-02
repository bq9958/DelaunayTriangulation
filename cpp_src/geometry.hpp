#pragma once

#include <cmath>
#include <array>

namespace dt {

struct Point {
    double x = 0.0;
    double y = 0.0;
};

inline double distance(const Point& a, const Point& b) {
    double dx = a.x - b.x, dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

inline double signed_area(const Point& p0, const Point& p1, const Point& p2) {
    return 0.5 * ((p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y));
}

inline double tri_area(const Point& p0, const Point& p1, const Point& p2) {
    return std::abs(signed_area(p0, p1, p2));
}

inline void barycentric(const Point& p0, const Point& p1, const Point& p2,
                        const Point& p, double& b0, double& b1, double& b2) {
    double area = signed_area(p0, p1, p2);
    b0 = signed_area(p, p1, p2) / area;
    b1 = signed_area(p0, p, p2) / area;
    b2 = signed_area(p0, p1, p) / area;
}

// Returns the local edge index (0,1,2) to walk toward, or -1 if p is inside the triangle.
// Edge i is opposite vertex i: edge 0 = (v1,v2), edge 1 = (v2,v0), edge 2 = (v0,v1).
inline int in_triangle(const Point& p0, const Point& p1, const Point& p2, const Point& p) {
    double b0, b1, b2;
    barycentric(p0, p1, p2, p, b0, b1, b2);

    int neg_count = 0;
    int neg_indices[3];
    double bary[3] = {b0, b1, b2};
    for (int i = 0; i < 3; i++) {
        if (bary[i] < 0.0) {
            neg_indices[neg_count++] = i;
        }
    }

    if (neg_count == 0) return -1;        // inside
    if (neg_count == 1) return neg_indices[0];
    // Two negative: pick one randomly (avoid deterministic cycling)
    return neg_indices[std::rand() % 2];
}

inline double circumcircle_radius(const Point& p0, const Point& p1, const Point& p2) {
    double l0 = distance(p0, p1);
    double l1 = distance(p1, p2);
    double l2 = distance(p2, p0);
    double area = tri_area(p0, p1, p2);
    return l0 * l1 * l2 / (4.0 * area);
}

inline Point circumcircle_center(const Point& p0, const Point& p1, const Point& p2) {
    double ax = p0.x, ay = p0.y;
    double bx = p1.x, by = p1.y;
    double cx = p2.x, cy = p2.y;

    double D = 2.0 * ((cx - ax) * (by - ay) - (bx - ax) * (cy - ay));
    double ux = ((by - ay) * (cx * cx - ax * ax + cy * cy - ay * ay)
               - (cy - ay) * (bx * bx - ax * ax + by * by - ay * ay)) / D;
    double uy = ((bx - ax) * (cx * cx - ax * ax + cy * cy - ay * ay)
               - (cx - ax) * (bx * bx - ax * ax + by * by - ay * ay)) / D;
    return {ux, uy};
}

inline bool in_circumcircle(const Point& p0, const Point& p1, const Point& p2, const Point& p) {
    // Determinant test: positive if p is inside circumcircle of CCW triangle (p0,p1,p2)
    double ax = p0.x - p.x, ay = p0.y - p.y;
    double bx = p1.x - p.x, by = p1.y - p.y;
    double cx = p2.x - p.x, cy = p2.y - p.y;

    double det = ax * (by * (cx * cx + cy * cy) - cy * (bx * bx + by * by))
               - ay * (bx * (cx * cx + cy * cy) - cx * (bx * bx + by * by))
               + (ax * ax + ay * ay) * (bx * cy - by * cx);

    return det > 0.0;
}

} // namespace dt
