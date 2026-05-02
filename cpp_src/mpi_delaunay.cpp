#include "mpi_delaunay.hpp"
#include <algorithm>
#include <numeric>
#include <iostream>
#include <cstdlib>

MPIDelaunay::MPIDelaunay(MPI_Comm comm) : comm_(comm) {
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);
}

void MPIDelaunay::triangulate_random(int num_points,
                                      double xmin, double xmax,
                                      double ymin, double ymax,
                                      unsigned seed) {
    std::vector<Point> all_points;
    if (rank_ == 0) {
        std::srand(seed ? seed : static_cast<unsigned>(num_points));
        all_points.resize(num_points);
        for (int i = 0; i < num_points; i++) {
            double rx = static_cast<double>(std::rand()) / RAND_MAX;
            double ry = static_cast<double>(std::rand()) / RAND_MAX;
            all_points[i] = {xmin + rx * (xmax - xmin), ymin + ry * (ymax - ymin)};
        }
    }
    triangulate(all_points, xmin, xmax, ymin, ymax);
}

MPIDelaunay::Subdomain MPIDelaunay::decompose(const std::vector<Point>& all_points,
                                               double xmin, double xmax,
                                               double ymin, double ymax) {
    Subdomain sub;
    int n = static_cast<int>(all_points.size());

    // Broadcast point count
    MPI_Bcast(&n, 1, MPI_INT, 0, comm_);

    // Sort points by x on rank 0 and broadcast
    std::vector<Point> sorted_pts(n);
    if (rank_ == 0) {
        std::vector<int> idx(n);
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b) {
            return all_points[a].x < all_points[b].x;
        });
        for (int i = 0; i < n; i++) sorted_pts[i] = all_points[idx[i]];
    }
    MPI_Bcast(sorted_pts.data(), n * 2, MPI_DOUBLE, 0, comm_);

    // Determine core region for this rank
    double domain_width = xmax - xmin;
    double slice = domain_width / size_;
    sub.xlo = xmin + rank_ * slice;
    sub.xhi = xmin + (rank_ + 1) * slice;
    sub.ylo = ymin;
    sub.yhi = ymax;

    // Guard band
    double guard = guard_fraction_ * domain_width;
    sub.xlo_ext = (rank_ == 0)        ? sub.xlo : sub.xlo - guard;
    sub.xhi_ext = (rank_ == size_ - 1) ? sub.xhi : sub.xhi + guard;

    // Collect points in extended region
    for (int i = 0; i < n; i++) {
        if (sorted_pts[i].x >= sub.xlo_ext && sorted_pts[i].x <= sub.xhi_ext) {
            sub.points.push_back(sorted_pts[i]);
        }
    }

    return sub;
}

void MPIDelaunay::merge(const Mesh& local, const Subdomain& sub) {
    // Count triangles whose centroid falls in this rank's core region
    std::vector<std::array<int,3>> core_tris;
    for (int i = 0; i < local.num_triangles; i++) {
        const auto& t = local.triangles[i];
        double cx = (local.vertices[t[0]].x + local.vertices[t[1]].x + local.vertices[t[2]].x) / 3.0;
        if (cx >= sub.xlo && cx < sub.xhi) {
            core_tris.push_back(t);
        }
    }

    // Gather triangle counts on rank 0
    int local_ntri = static_cast<int>(core_tris.size());
    int local_nver = static_cast<int>(local.vertices.size());
    std::vector<int> all_ntri(size_), all_nver(size_);
    MPI_Gather(&local_ntri, 1, MPI_INT, all_ntri.data(), 1, MPI_INT, 0, comm_);
    MPI_Gather(&local_nver, 1, MPI_INT, all_nver.data(), 1, MPI_INT, 0, comm_);

    // Gather all vertices (as flat doubles)
    std::vector<double> local_verts_flat(local_nver * 2);
    for (int i = 0; i < local_nver; i++) {
        local_verts_flat[2*i]   = local.vertices[i].x;
        local_verts_flat[2*i+1] = local.vertices[i].y;
    }

    // Gather all core triangles (as flat ints)
    std::vector<int> local_tris_flat(local_ntri * 3);
    for (int i = 0; i < local_ntri; i++) {
        local_tris_flat[3*i]   = core_tris[i][0];
        local_tris_flat[3*i+1] = core_tris[i][1];
        local_tris_flat[3*i+2] = core_tris[i][2];
    }

    if (rank_ == 0) {
        // Compute displacements
        std::vector<int> vdispl(size_), tdispl(size_);
        std::vector<int> vcount(size_), tcount(size_);
        int total_verts = 0, total_tris = 0;
        for (int r = 0; r < size_; r++) {
            vcount[r] = all_nver[r] * 2;
            tcount[r] = all_ntri[r] * 3;
            vdispl[r] = total_verts * 2;
            tdispl[r] = total_tris * 3;
            total_verts += all_nver[r];
            total_tris += all_ntri[r];
        }

        std::vector<double> all_verts(total_verts * 2);
        std::vector<int>    all_tris(total_tris * 3);

        MPI_Gatherv(local_verts_flat.data(), local_nver * 2, MPI_DOUBLE,
                    all_verts.data(), vcount.data(), vdispl.data(), MPI_DOUBLE, 0, comm_);
        MPI_Gatherv(local_tris_flat.data(), local_ntri * 3, MPI_INT,
                    all_tris.data(), tcount.data(), tdispl.data(), MPI_INT, 0, comm_);

        // Build global mesh by concatenating (with vertex index offsets)
        global_mesh_.vertices.resize(total_verts);
        for (int i = 0; i < total_verts; i++) {
            global_mesh_.vertices[i] = {all_verts[2*i], all_verts[2*i+1]};
        }

        global_mesh_.triangles.resize(total_tris);
        global_mesh_.neighbors.resize(total_tris, {-1, -1, -1});
        int vert_offset = 0;
        int tri_idx = 0;
        for (int r = 0; r < size_; r++) {
            for (int i = 0; i < all_ntri[r]; i++) {
                int base = tdispl[r] / 3 + i;
                global_mesh_.triangles[tri_idx] = {
                    all_tris[3*base]   + vert_offset,
                    all_tris[3*base+1] + vert_offset,
                    all_tris[3*base+2] + vert_offset
                };
                tri_idx++;
            }
            vert_offset += all_nver[r];
        }
        global_mesh_.num_triangles = total_tris;

        std::cout << "[MPI] Merged mesh: " << total_verts << " vertices, "
                  << total_tris << " triangles\n";
    } else {
        MPI_Gatherv(local_verts_flat.data(), local_nver * 2, MPI_DOUBLE,
                    nullptr, nullptr, nullptr, MPI_DOUBLE, 0, comm_);
        MPI_Gatherv(local_tris_flat.data(), local_ntri * 3, MPI_INT,
                    nullptr, nullptr, nullptr, MPI_INT, 0, comm_);
    }
}

void MPIDelaunay::triangulate(const std::vector<Point>& all_points,
                               double xmin, double xmax,
                               double ymin, double ymax) {
    auto sub = decompose(all_points, xmin, xmax, ymin, ymax);

    if (rank_ == 0) {
        std::cout << "[MPI] Rank " << rank_ << ": " << sub.points.size()
                  << " points (core [" << sub.xlo << ", " << sub.xhi << "])\n";
    }

    // Local triangulation
    Delaunay local_dt;
    local_dt.triangulate(sub.points);

    if (rank_ == 0) {
        std::cout << "[MPI] Rank " << rank_ << ": "
                  << local_dt.get_mesh().num_triangles << " local triangles\n";
    }

    merge(local_dt.get_mesh(), sub);
}
