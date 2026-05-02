#include "delaunay.hpp"
#include "mesh_io.hpp"
#include <iostream>
#include <chrono>

int main(int argc, char* argv[]) {
    int num_points = 200000;
    if (argc > 1) num_points = std::atoi(argv[1]);

    std::cout << "Delaunay triangulation of " << num_points << " random points\n";

    Delaunay dt;
    auto t0 = std::chrono::high_resolution_clock::now();
    dt.triangulate_random(num_points, 0, 1, 0, 1);
    auto t1 = std::chrono::high_resolution_clock::now();

    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "Time: " << elapsed << " s\n";
    std::cout << "Vertices: " << dt.get_mesh().vertices.size()
              << "  Triangles: " << dt.get_mesh().num_triangles << "\n";

    MeshIO::write_triangulation(dt.get_mesh(), "triangulation.mesh");
    MeshIO::write_vtk(dt.get_mesh(), "triangulation.vtk");
    return 0;
}
