#pragma once

#include "mesh.hpp"
#include <string>

// Pure-STL replacement for libmesh6. Handles ASCII .mesh format (1-indexed on disk,
// 0-indexed in memory).
namespace MeshIO {
    // Read a .mesh file. Returns empty mesh on failure.
    Mesh read(const std::string& path);

    // Write mesh to a .mesh file.
    bool write(const Mesh& mesh, const std::string& path);

    // Write mesh in Vizir-compatible format (includes all vertices and triangles).
    bool write_triangulation(const Mesh& mesh, const std::string& path, int num_frontier_pts = 0);

    // Write mesh as VTK Legacy ASCII (readable by ParaView).
    bool write_vtk(const Mesh& mesh, const std::string& path);
}
