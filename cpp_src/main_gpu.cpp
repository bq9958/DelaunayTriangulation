#include "delaunay.hpp"
#include "metal_backend.hpp"
#include "mesh_io.hpp"
#include <iostream>
#include <chrono>
#include <vector>

// Demonstrates GPU-accelerated batch circumcircle testing.
// Builds a serial Delaunay triangulation, then verifies the Delaunay property
// using the Metal GPU backend for batch circumcircle checks.
int main(int argc, char* argv[]) {
    int num_points = 200000;
    if (argc > 1) num_points = std::atoi(argv[1]);

    std::cout << "Delaunay triangulation with Metal GPU verification\n";
    std::cout << "Points: " << num_points << "\n";

    // Step 1: Build triangulation (serial)
    Delaunay dt;
    auto t0 = std::chrono::high_resolution_clock::now();
    dt.triangulate_random(num_points, 0, 1, 0, 1);
    auto t1 = std::chrono::high_resolution_clock::now();

    double build_time = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "Build time: " << build_time << " s\n";

    const Mesh& mesh = dt.get_mesh();
    std::cout << "Vertices: " << mesh.vertices.size()
              << "  Triangles: " << mesh.num_triangles << "\n";

    // Step 2: Verify Delaunay property using GPU batch circumcircle test
    MetalBackend metal;

    // For each vertex, check that it does not lie inside any non-adjacent circumcircle.
    // We batch all triangles against each test point.
    std::vector<std::array<int,3>> all_tris(mesh.triangles.begin(),
                                             mesh.triangles.begin() + mesh.num_triangles);

    int violations = 0;
    auto t2 = std::chrono::high_resolution_clock::now();

    // Sample a subset of vertices for verification (checking all is O(V*T))
    int sample_count = std::min(1000, static_cast<int>(mesh.vertices.size()));
    for (int vi = 0; vi < sample_count; vi++) {
        int vert = vi * static_cast<int>(mesh.vertices.size()) / sample_count;
        auto results = metal.batch_in_circumcircle(mesh.vertices, all_tris,
                                                    mesh.vertices[vert]);
        for (int ti = 0; ti < mesh.num_triangles; ti++) {
            if (!results[ti]) continue;
            // Check that vert is not a vertex of this triangle
            const auto& t = all_tris[ti];
            if (t[0] != vert && t[1] != vert && t[2] != vert) {
                violations++;
            }
        }
    }

    auto t3 = std::chrono::high_resolution_clock::now();
    double verify_time = std::chrono::duration<double>(t3 - t2).count();

    std::cout << "Verification time (" << sample_count << " samples): "
              << verify_time << " s"
              << (metal.is_available() ? " [GPU]" : " [CPU fallback]") << "\n";
    std::cout << "Delaunay violations: " << violations << "\n";

    MeshIO::write_triangulation(mesh, "triangulation_gpu.mesh");
    MeshIO::write_vtk(mesh, "triangulation_gpu.vtk");
    return 0;
}
