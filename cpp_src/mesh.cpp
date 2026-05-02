#include "mesh.hpp"
#include <cmath>
#include <algorithm>

using namespace dt;

void Mesh::compute_bounding_box() {
    box = {1e30, -1e30, 1e30, -1e30};
    for (size_t i = 4; i < vertices.size(); i++) {
        box[0] = std::min(box[0], vertices[i].x);
        box[1] = std::max(box[1], vertices[i].x);
        box[2] = std::min(box[2], vertices[i].y);
        box[3] = std::max(box[3], vertices[i].y);
    }
    double dx = box[1] - box[0];
    double dy = box[3] - box[2];
    box[0] -= 0.5 * dx; box[1] += 0.5 * dx;
    box[2] -= 0.5 * dy; box[3] += 0.5 * dy;

    // Set bounding box corners (top-left, bottom-left, bottom-right, top-right)
    vertices[0] = {box[0], box[3]};
    vertices[1] = {box[0], box[2]};
    vertices[2] = {box[1], box[2]};
    vertices[3] = {box[1], box[3]};
}

void Mesh::quality(std::vector<double>& qual, int mode) const {
    const double alpha1 = 1.0 / (4.0 * std::sqrt(3.0));
    const double alpha2 = std::sqrt(3.0) / 6.0;

    qual.resize(num_triangles);
    for (int i = 0; i < num_triangles; i++) {
        const auto& t = triangles[i];
        const Point& p0 = vertices[t[0]];
        const Point& p1 = vertices[t[1]];
        const Point& p2 = vertices[t[2]];

        double l0 = distance(p0, p1);
        double l1 = distance(p1, p2);
        double l2 = distance(p2, p0);
        double area = tri_area(p0, p1, p2);

        if (mode == 1) {
            qual[i] = alpha1 * (l0 * l0 + l1 * l1 + l2 * l2) / area;
        } else {
            double rho = 2.0 * area / (l0 + l1 + l2);
            qual[i] = alpha2 * std::max({l0, l1, l2}) / rho;
        }
    }
}

void Mesh::reserve(int num_verts, int num_tris) {
    vertices.reserve(4 + num_verts);
    triangles.reserve(num_tris);
    neighbors.reserve(num_tris);
}
