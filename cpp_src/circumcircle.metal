#include <metal_stdlib>
using namespace metal;

// Batch circumcircle test on GPU.
// For each candidate triangle, determine whether query_point lies inside its circumcircle.
// Uses the determinant formulation (positive ⟺ inside for CCW triangles).
kernel void batch_circumcircle_test(
    device const float2* vertices        [[buffer(0)]],
    device const int*    tri_indices     [[buffer(1)]],  // packed [v0,v1,v2, v0,v1,v2, ...]
    device const float2& query_point    [[buffer(2)]],
    device int*          results         [[buffer(3)]],  // 1 = inside, 0 = outside
    uint tid [[thread_position_in_grid]])
{
    int i0 = tri_indices[tid * 3];
    int i1 = tri_indices[tid * 3 + 1];
    int i2 = tri_indices[tid * 3 + 2];

    float2 a = vertices[i0];
    float2 b = vertices[i1];
    float2 c = vertices[i2];
    float2 d = query_point;

    float ax = a.x - d.x, ay = a.y - d.y;
    float bx = b.x - d.x, by = b.y - d.y;
    float cx = c.x - d.x, cy = c.y - d.y;

    float a2 = ax * ax + ay * ay;
    float b2 = bx * bx + by * by;
    float c2 = cx * cx + cy * cy;

    float det = ax * (by * c2 - cy * b2)
              - ay * (bx * c2 - cx * b2)
              + a2 * (bx * cy - by * cx);

    results[tid] = (det > 0.0f) ? 1 : 0;
}
