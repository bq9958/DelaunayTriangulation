#include "mpi_delaunay.hpp"
#include "mesh_io.hpp"
#include <iostream>
#include <chrono>
#include <mpi.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int num_points = 200000;
    if (argc > 1) num_points = std::atoi(argv[1]);

    MPIDelaunay mpi_dt;

    if (mpi_dt.rank() == 0)
        std::cout << "MPI Delaunay triangulation of " << num_points << " points\n";

    auto t0 = std::chrono::high_resolution_clock::now();
    mpi_dt.triangulate_random(num_points, 0, 1, 0, 1);
    auto t1 = std::chrono::high_resolution_clock::now();

    if (mpi_dt.rank() == 0) {
        double elapsed = std::chrono::duration<double>(t1 - t0).count();
        std::cout << "Time: " << elapsed << " s\n";
        MeshIO::write_triangulation(mpi_dt.get_mesh(), "triangulation_mpi.mesh");
        MeshIO::write_vtk(mpi_dt.get_mesh(), "triangulation_mpi.vtk");
    }

    MPI_Finalize();
    return 0;
}
