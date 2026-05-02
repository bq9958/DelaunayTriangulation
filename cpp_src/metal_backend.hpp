#pragma once

#include "geometry.hpp"
#include <vector>
#include <array>
#include <memory>

// Metal GPU backend for batch circumcircle tests.
// Falls back to CPU if Metal is unavailable.
class MetalBackend {
public:
    MetalBackend();
    ~MetalBackend();

    bool is_available() const;

    // Test whether query_point is inside the circumcircle of each candidate triangle.
    // vertices: full vertex array. candidates: list of {v0,v1,v2} index triples.
    // Returns a vector of bools (true = inside circumcircle).
    std::vector<bool> batch_in_circumcircle(
        const std::vector<dt::Point>& vertices,
        const std::vector<std::array<int,3>>& candidates,
        const dt::Point& query_point);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
