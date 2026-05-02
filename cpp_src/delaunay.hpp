#pragma once

#include "mesh.hpp"
#include "edge_map.hpp"
#include <vector>
#include <functional>

class Delaunay {
public:
    // Build Delaunay triangulation of num_points random points in [xmin,xmax] x [ymin,ymax].
    void triangulate_random(int num_points,
                            double xmin = 0, double xmax = 1,
                            double ymin = 0, double ymax = 1,
                            unsigned seed = 0);

    // Build Delaunay triangulation from a pre-loaded set of points.
    // points should NOT include the 4 bounding-box vertices; they will be prepended.
    void triangulate(const std::vector<Point>& points);

    Mesh& get_mesh() { return mesh_; }
    const Mesh& get_mesh() const { return mesh_; }

private:
    Mesh mesh_;
    EdgeMap edge_map_;
    int last_located_tri_ = 0;

    // Initialize mesh with bounding box and 2 starter triangles.
    void init_bounding_box(int num_points);

    // Insert point at vertex index vert_id (must already be in mesh_.vertices).
    void insert_point(int vert_id);

    // Walk from tri toward the triangle containing p. Returns triangle index.
    int locate(int start_tri, const Point& p) const;

    // Find the two vertices shared between triangles i and j.
    // Returns (global_v1, global_v2, local_in_j_1, local_in_j_2).
    static std::array<int,4> common_edge(const Mesh& m, int i, int j);
};
