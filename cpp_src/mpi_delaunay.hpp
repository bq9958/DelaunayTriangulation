#pragma once

#include "delaunay.hpp"
#include <mpi.h>
#include <vector>

// MPI-parallel Delaunay triangulation via recursive coordinate bisection.
//
// Strategy:
//   1. Rank 0 generates (or reads) the full point set.
//   2. Points are sorted by x-coordinate and split into P contiguous chunks.
//   3. Each chunk is expanded by a guard band (overlap) so that boundary
//      triangles remain valid.
//   4. Each rank builds a local Delaunay triangulation on its chunk.
//   5. The global mesh is gathered on rank 0 by merging local results,
//      keeping only triangles whose centroid lies in the rank's core region.
class MPIDelaunay {
public:
    explicit MPIDelaunay(MPI_Comm comm = MPI_COMM_WORLD);

    // Triangulate num_points random points. Result available on rank 0.
    void triangulate_random(int num_points,
                            double xmin = 0, double xmax = 1,
                            double ymin = 0, double ymax = 1,
                            unsigned seed = 0);

    // Triangulate a given point set (only rank 0 needs to provide it).
    void triangulate(const std::vector<Point>& all_points,
                     double xmin, double xmax, double ymin, double ymax);

    Mesh& get_mesh() { return global_mesh_; }
    int rank() const { return rank_; }

private:
    MPI_Comm comm_;
    int rank_, size_;
    Mesh global_mesh_;

    // Fraction of subdomain width used as overlap on each side.
    static constexpr double guard_fraction_ = 0.05;

    struct Subdomain {
        double xlo, xhi;           // core region
        double xlo_ext, xhi_ext;   // extended region (with guard band)
        double ylo, yhi;
        std::vector<Point> points; // points in extended region
    };

    Subdomain decompose(const std::vector<Point>& all_points,
                        double xmin, double xmax,
                        double ymin, double ymax);

    void merge(const Mesh& local, const Subdomain& sub);
};
