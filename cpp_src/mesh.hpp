#pragma once

#include "geometry.hpp"
#include <vector>
#include <array>

using dt::Point;

class Mesh {
public:
    int dim = 2;

    std::vector<Point>             vertices;    // 0-indexed
    std::vector<std::array<int,3>> triangles;   // vertex indices (0-indexed)
    std::vector<std::array<int,3>> neighbors;   // neighbor triangle indices, -1 = none
    std::vector<std::array<int,2>> boundary_edges;
    std::array<double,4>           box = {1e30, -1e30, 1e30, -1e30}; // xmin,xmax,ymin,ymax

    int num_triangles = 0; // active triangle count (may differ from triangles.size())

    // Compute bounding box from vertices[4..] and set vertices[0..3] to enlarged corners.
    void compute_bounding_box();

    // Compute triangle quality. mode 1: circumradius/area, mode 2: longest_edge/inradius.
    void quality(std::vector<double>& qual, int mode) const;

    // Reserve storage for expected sizes.
    void reserve(int num_verts, int num_tris);
};
