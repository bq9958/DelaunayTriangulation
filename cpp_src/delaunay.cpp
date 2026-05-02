#include "delaunay.hpp"
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <cassert>

using namespace dt;

// ---------- helpers ----------

static std::vector<Point> generate_random_points(int n, double xmin, double xmax,
                                                  double ymin, double ymax, unsigned seed) {
    std::srand(seed ? seed : static_cast<unsigned>(n));
    std::vector<Point> pts(n);
    for (int i = 0; i < n; i++) {
        double rx = static_cast<double>(std::rand()) / RAND_MAX;
        double ry = static_cast<double>(std::rand()) / RAND_MAX;
        pts[i] = {xmin + rx * (xmax - xmin), ymin + ry * (ymax - ymin)};
    }
    return pts;
}

// Remove (v1,v2) or (v2,v1) from edge list.
static void remove_boundary_edge(std::vector<std::pair<int,int>>& edges, int v1, int v2) {
    edges.erase(std::remove_if(edges.begin(), edges.end(),
        [v1, v2](const std::pair<int,int>& e) {
            return (e.first == v1 && e.second == v2)
                || (e.first == v2 && e.second == v1);
        }), edges.end());
}

// Sort directed edges into a circular polygon.
static void sort_edges_circular(std::vector<std::pair<int,int>>& edges) {
    if (edges.size() <= 1) return;
    std::vector<std::pair<int,int>> sorted;
    sorted.reserve(edges.size());
    std::vector<bool> used(edges.size(), false);

    sorted.push_back(edges[0]);
    used[0] = true;
    for (size_t i = 1; i < edges.size(); i++) {
        int tail = sorted.back().second;
        for (size_t j = 1; j < edges.size(); j++) {
            if (!used[j] && edges[j].first == tail) {
                sorted.push_back(edges[j]);
                used[j] = true;
                break;
            }
        }
    }
    edges = std::move(sorted);
}

// ---------- Delaunay ----------

void Delaunay::triangulate_random(int num_points, double xmin, double xmax,
                                   double ymin, double ymax, unsigned seed) {
    auto pts = generate_random_points(num_points, xmin, xmax, ymin, ymax, seed);
    triangulate(pts);
}

void Delaunay::triangulate(const std::vector<Point>& points) {
    int n = static_cast<int>(points.size());

    // Prepare vertices: [0..3] = bounding box (filled later), [4..4+n-1] = input points
    mesh_.vertices.resize(4 + n);
    for (int i = 0; i < n; i++) {
        mesh_.vertices[4 + i] = points[i];
    }

    init_bounding_box(n);

    last_located_tri_ = 0;
    for (int i = 0; i < n; i++) {
        insert_point(4 + i);
    }
}

void Delaunay::init_bounding_box(int num_points) {
    int cap = 4 * num_points;

    mesh_.triangles.resize(cap, {0, 0, 0});
    mesh_.neighbors.resize(cap, {-1, -1, -1});
    mesh_.num_triangles = 0;

    mesh_.compute_bounding_box();

    //  0(TL)-----3(TR)
    //  |  tri0  / |
    //  |      /   |
    //  |    /     |
    //  |  / tri1  |
    //  1(BL)-----2(BR)
    mesh_.triangles[0] = {0, 1, 3};
    mesh_.triangles[1] = {1, 2, 3};
    mesh_.num_triangles = 2;

    // Edge 0 of tri0 = (v1,v2) = (1,3); edge 1 of tri1 = (v2,v0) = (3,1). They share edge (1,3).
    mesh_.neighbors[0] = { 1, -1, -1};
    mesh_.neighbors[1] = {-1,  0, -1};

    edge_map_.clear();
    // Store initial edges for neighbor lookup during insertion
    // tri0 edges: (1,3) edge0, (3,0) edge1, (0,1) edge2
    edge_map_.add(1, 3, 0, 0);
    edge_map_.add(3, 0, 0, 1);
    edge_map_.add(0, 1, 0, 2);
    // tri1 edges: (2,3) edge0, (3,1) edge1, (1,2) edge2
    edge_map_.add(2, 3, 1, 0);
    edge_map_.add(3, 1, 1, 1);
    edge_map_.add(1, 2, 1, 2);
}

int Delaunay::locate(int start_tri, const Point& p) const {
    int tri = start_tri;
    int max_steps = 2 * mesh_.num_triangles;
    for (int step = 0; step < max_steps; step++) {
        const auto& t = mesh_.triangles[tri];
        int dir = in_triangle(mesh_.vertices[t[0]], mesh_.vertices[t[1]],
                              mesh_.vertices[t[2]], p);
        if (dir == -1) return tri;  // found

        int next = mesh_.neighbors[tri][dir];
        if (next < 0) {
            std::cerr << "[ERROR] Location walk reached boundary\n";
            return tri;
        }
        tri = next;
    }
    std::cerr << "[ERROR] Location walk exceeded step limit\n";
    return tri;
}

std::array<int,4> Delaunay::common_edge(const Mesh& m, int i, int j) {
    // Returns {global_v1, global_v2, local_j_1, local_j_2}
    std::array<int,4> result = {-1, -1, -1, -1};
    int count = 0;
    for (int li = 0; li < 3 && count < 2; li++) {
        for (int lj = 0; lj < 3; lj++) {
            if (m.triangles[j][lj] == m.triangles[i][li]) {
                if (count == 0) {
                    result[0] = m.triangles[j][lj];
                    result[2] = lj;
                } else {
                    result[1] = m.triangles[j][lj];
                    result[3] = lj;
                }
                count++;
                break;
            }
        }
    }
    return result;
}

void Delaunay::insert_point(int vert_id) {
    const Point& pt = mesh_.vertices[vert_id];

    // --- 1. Locate containing triangle ---
    if (last_located_tri_ <= 0 || last_located_tri_ >= mesh_.num_triangles)
        last_located_tri_ = std::rand() % mesh_.num_triangles;
    int tri_loc = locate(last_located_tri_, pt);
    last_located_tri_ = tri_loc;

    // --- 2. Build cavity ---
    std::vector<int> cavity = {tri_loc};            // triangle indices in cavity
    std::vector<std::pair<int,int>> boundary;       // directed boundary edges

    const auto& t0 = mesh_.triangles[tri_loc];
    boundary.push_back({t0[1], t0[2]});  // edge 0
    boundary.push_back({t0[2], t0[0]});  // edge 1
    boundary.push_back({t0[0], t0[1]});  // edge 2

    int ptr_cur = 0;
    while (ptr_cur < static_cast<int>(cavity.size())) {
        int depile = cavity[ptr_cur++];
        for (int edge_local = 0; edge_local < 3; edge_local++) {
            int neighbor = mesh_.neighbors[depile][edge_local];
            if (neighbor < 0) continue;
            // Check if already in cavity
            if (std::find(cavity.begin(), cavity.end(), neighbor) != cavity.end())
                continue;

            const auto& tn = mesh_.triangles[neighbor];
            if (in_circumcircle(mesh_.vertices[tn[0]], mesh_.vertices[tn[1]],
                                mesh_.vertices[tn[2]], pt)) {
                cavity.push_back(neighbor);

                // Find common edge between depile and neighbor
                auto ce = common_edge(mesh_, depile, neighbor);
                int cv1 = ce[0], cv2 = ce[1];
                int cl1 = ce[2], cl2 = ce[3];

                // Remove common edge from boundary
                remove_boundary_edge(boundary, cv1, cv2);

                // Add the other two edges of neighbor
                int reste = 3 - cl1 - cl2;
                int vmin_l = std::min(cl1, cl2), vmax_l = std::max(cl1, cl2);
                int vmin_g = mesh_.triangles[neighbor][vmin_l];
                int vmax_g = mesh_.triangles[neighbor][vmax_l];
                int vreste_g = mesh_.triangles[neighbor][reste];

                // Add in correct winding order
                if (reste == 1) {
                    boundary.push_back({vmin_g, vreste_g});
                    boundary.push_back({vreste_g, vmax_g});
                } else { // reste == 0 or 2
                    boundary.push_back({vmax_g, vreste_g});
                    boundary.push_back({vreste_g, vmin_g});
                }
            }
        }
    }

    // --- 3. Sort boundary edges into a polygon ---
    sort_edges_circular(boundary);

    // --- 4. Re-triangulate cavity ---
    int num_new = static_cast<int>(boundary.size());
    int num_old = static_cast<int>(cavity.size());

    // Slots for new triangles: reuse cavity slots first, then append
    std::vector<int> new_tri_ids;
    new_tri_ids.reserve(num_new);
    for (int i = 0; i < num_new; i++) {
        if (i < num_old) {
            new_tri_ids.push_back(cavity[i]);
        } else {
            int slot = mesh_.num_triangles + (i - num_old);
            new_tri_ids.push_back(slot);
        }
    }
    mesh_.num_triangles += (num_new - num_old);

    // Create new triangles and update edge map
    for (int i = 0; i < num_new; i++) {
        int slot = new_tri_ids[i];
        int ev1 = boundary[i].first;
        int ev2 = boundary[i].second;

        mesh_.triangles[slot] = {vert_id, ev1, ev2};

        // Register the two radial edges in the edge map
        // Edge 2: (v0, v1) = (vert_id, ev1)
        edge_map_.add(vert_id, ev1, slot, 2);
        // Edge 1: (v2, v0) = (ev2, vert_id)
        edge_map_.add(ev2, vert_id, slot, 1);

        // Update edge map for the boundary edge (ev1, ev2) = edge 0
        auto* info = edge_map_.find_directed(ev1, ev2);
        if (info) {
            info->tri = slot;
            info->edge = 0;
        }

        // Find external neighbor across boundary edge (reverse direction)
        auto* ext = edge_map_.find_directed(ev2, ev1);
        int ext_neighbor = ext ? ext->tri : -1;

        if (ext_neighbor >= 0) {
            // Find which local edge of ext_neighbor shares (ev2, ev1)
            auto ce = common_edge(mesh_, slot, ext_neighbor);
            int local_in_ext = 3 - ce[2] - ce[3];
            mesh_.neighbors[ext_neighbor][local_in_ext] = slot;
            mesh_.neighbors[slot][0] = ext_neighbor;
        } else {
            mesh_.neighbors[slot][0] = -1;
        }
    }

    // Connect adjacent cavity triangles (radial neighbors around the new point)
    for (int i = 0; i < num_new; i++) {
        int prev = (i == 0) ? num_new - 1 : i - 1;
        int next = (i == num_new - 1) ? 0 : i + 1;
        mesh_.neighbors[new_tri_ids[i]][2] = new_tri_ids[prev];  // edge 2 neighbor
        mesh_.neighbors[new_tri_ids[i]][1] = new_tri_ids[next];  // edge 1 neighbor
    }
}
