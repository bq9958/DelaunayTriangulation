#include "mesh_io.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cctype>

static std::string trim(const std::string& s) {
    auto start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

// Skip blank lines and comments, return next non-empty trimmed line.
static bool next_line(std::ifstream& f, std::string& line) {
    while (std::getline(f, line)) {
        line = trim(line);
        if (!line.empty() && line[0] != '#') return true;
    }
    return false;
}

Mesh MeshIO::read(const std::string& path) {
    Mesh mesh;
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "  ## ERROR: Cannot open " << path << "\n";
        return mesh;
    }

    std::string line;
    while (next_line(file, line)) {
        if (line.find("MeshVersionFormatted") != std::string::npos) {
            continue;
        }
        else if (line.find("Dimension") != std::string::npos) {
            // Dimension might be on this line or the next
            std::istringstream iss(line);
            std::string kw;
            iss >> kw;
            if (!(iss >> mesh.dim)) {
                next_line(file, line);
                mesh.dim = std::stoi(line);
            }
        }
        else if (line == "Vertices") {
            next_line(file, line);
            int n = std::stoi(line);
            mesh.vertices.resize(n);
            for (int i = 0; i < n; i++) {
                next_line(file, line);
                std::istringstream iss(line);
                int ref;
                iss >> mesh.vertices[i].x >> mesh.vertices[i].y >> ref;
            }
        }
        else if (line == "Triangles") {
            next_line(file, line);
            int n = std::stoi(line);
            mesh.triangles.resize(n);
            mesh.neighbors.resize(n, {-1, -1, -1});
            mesh.num_triangles = n;
            for (int i = 0; i < n; i++) {
                next_line(file, line);
                std::istringstream iss(line);
                int ref;
                iss >> mesh.triangles[i][0] >> mesh.triangles[i][1]
                    >> mesh.triangles[i][2] >> ref;
                // Convert 1-indexed → 0-indexed
                mesh.triangles[i][0]--;
                mesh.triangles[i][1]--;
                mesh.triangles[i][2]--;
            }
        }
        else if (line == "Edges") {
            next_line(file, line);
            int n = std::stoi(line);
            mesh.boundary_edges.resize(n);
            for (int i = 0; i < n; i++) {
                next_line(file, line);
                std::istringstream iss(line);
                int ref;
                iss >> mesh.boundary_edges[i][0] >> mesh.boundary_edges[i][1] >> ref;
                mesh.boundary_edges[i][0]--;
                mesh.boundary_edges[i][1]--;
            }
        }
        else if (line == "End") {
            break;
        }
    }

    std::cout << "  File " << path << " opened  Dimension " << mesh.dim
              << "  Vertices " << mesh.vertices.size()
              << "  Triangles " << mesh.num_triangles << "\n";
    return mesh;
}

bool MeshIO::write(const Mesh& mesh, const std::string& path) {
    std::ofstream file(path);
    if (!file.is_open()) {
        std::cerr << "  ## ERROR: Cannot create " << path << "\n";
        return false;
    }

    file << "MeshVersionFormatted 2\n\nDimension 2\n\n";

    file << "Vertices\n" << mesh.vertices.size() << "\n";
    for (const auto& v : mesh.vertices) {
        file << v.x << " " << v.y << " 0\n";  // 1-indexed not needed for coords
    }

    file << "\nTriangles\n" << mesh.num_triangles << "\n";
    for (int i = 0; i < mesh.num_triangles; i++) {
        const auto& t = mesh.triangles[i];
        file << (t[0]+1) << " " << (t[1]+1) << " " << (t[2]+1) << " 1\n";
    }

    if (!mesh.boundary_edges.empty()) {
        file << "\nEdges\n" << mesh.boundary_edges.size() << "\n";
        for (const auto& e : mesh.boundary_edges) {
            file << (e[0]+1) << " " << (e[1]+1) << " 1\n";
        }
    }

    file << "\nEnd\n";
    return true;
}

bool MeshIO::write_triangulation(const Mesh& mesh, const std::string& path, int num_frontier_pts) {
    std::ofstream file(path);
    if (!file.is_open()) {
        std::cerr << "  ## ERROR: Cannot create " << path << "\n";
        return false;
    }

    int total_verts = static_cast<int>(mesh.vertices.size());

    file << "MeshVersionFormatted 2\n\nDimension 2\n\n";

    file << "Vertices\n" << total_verts << "\n";
    for (int i = 0; i < total_verts; i++) {
        file << mesh.vertices[i].x << " " << mesh.vertices[i].y << " 0\n";
    }

    file << "\n\nTriangles\n" << mesh.num_triangles << "\n";
    for (int i = 0; i < mesh.num_triangles; i++) {
        const auto& t = mesh.triangles[i];
        file << (t[0]+1) << " " << (t[1]+1) << " " << (t[2]+1) << " 1\n";
    }

    if (num_frontier_pts > 0) {
        int num_edges = num_frontier_pts + 4;
        file << "\n\nEdges\n" << num_edges << "\n";
        // Bounding box edges (0-indexed vertices 0,1,2,3 → 1-indexed 1,2,3,4)
        file << "1 2 1\n2 3 1\n3 4 1\n4 1 1\n";
        for (int i = 1; i <= num_frontier_pts - 1; i++) {
            file << (4 + i) << " " << (5 + i) << " 2\n";
        }
        file << (4 + num_frontier_pts) << " 5 2\n";
    }

    file << "End\n";
    std::cout << "[Output] Triangulation written to " << path << "\n";
    return true;
}

bool MeshIO::write_vtk(const Mesh& mesh, const std::string& path) {
    std::ofstream file(path);
    if (!file.is_open()) {
        std::cerr << "  ## ERROR: Cannot create " << path << "\n";
        return false;
    }

    int nv = static_cast<int>(mesh.vertices.size());
    int nt = mesh.num_triangles;

    // Header
    file << "# vtk DataFile Version 3.0\n";
    file << "Delaunay Triangulation\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    // Points (2D → z = 0)
    file << "POINTS " << nv << " double\n";
    for (int i = 0; i < nv; i++) {
        file << mesh.vertices[i].x << " " << mesh.vertices[i].y << " 0.0\n";
    }

    // Cells (triangles, each line: 3 v0 v1 v2)
    file << "CELLS " << nt << " " << nt * 4 << "\n";
    for (int i = 0; i < nt; i++) {
        const auto& t = mesh.triangles[i];
        file << "3 " << t[0] << " " << t[1] << " " << t[2] << "\n";
    }

    // Cell types (5 = VTK_TRIANGLE)
    file << "CELL_TYPES " << nt << "\n";
    for (int i = 0; i < nt; i++) {
        file << "5\n";
    }

    std::cout << "[Output] VTK written to " << path << "\n";
    return true;
}
