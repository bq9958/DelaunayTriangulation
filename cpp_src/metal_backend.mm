#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#include "metal_backend.hpp"
#include <iostream>
#include <vector>

struct MetalBackend::Impl {
    id<MTLDevice>              device = nil;
    id<MTLCommandQueue>        queue  = nil;
    id<MTLComputePipelineState> pipeline = nil;
    bool available = false;

    Impl() {
        @autoreleasepool {
            device = MTLCreateSystemDefaultDevice();
            if (!device) {
                std::cerr << "[Metal] No GPU device found\n";
                return;
            }
            queue = [device newCommandQueue];

            // Load the Metal shader library from the .metallib compiled alongside this binary,
            // or compile from source at runtime.
            NSError* error = nil;

            // Try loading a pre-compiled .metallib next to the executable
            NSString* libPath = [[NSBundle mainBundle] pathForResource:@"circumcircle"
                                                               ofType:@"metallib"];
            id<MTLLibrary> library = nil;
            if (libPath) {
                NSURL* libURL = [NSURL fileURLWithPath:libPath];
                library = [device newLibraryWithURL:libURL error:&error];
            }

            // Fallback: compile from source at runtime (works without full Xcode)
            if (!library) {
                NSString* source = nil;
                // Use the compile-time path injected by CMake
#ifdef METAL_SHADER_PATH
                source = [NSString stringWithContentsOfFile:@METAL_SHADER_PATH
                                                  encoding:NSUTF8StringEncoding
                                                     error:&error];
#endif
                // Also try common relative paths
                if (!source) {
                    for (NSString* p in @[@"circumcircle.metal",
                                          @"cpp_src/circumcircle.metal"]) {
                        source = [NSString stringWithContentsOfFile:p
                                                          encoding:NSUTF8StringEncoding
                                                             error:&error];
                        if (source) break;
                    }
                }
                if (source) {
                    library = [device newLibraryWithSource:source options:nil error:&error];
                }
            }

            if (!library) {
                std::cerr << "[Metal] Cannot load shader library: "
                          << [[error localizedDescription] UTF8String] << "\n";
                return;
            }

            id<MTLFunction> func = [library newFunctionWithName:@"batch_circumcircle_test"];
            if (!func) {
                std::cerr << "[Metal] Kernel function not found\n";
                return;
            }

            pipeline = [device newComputePipelineStateWithFunction:func error:&error];
            if (!pipeline) {
                std::cerr << "[Metal] Pipeline creation failed: "
                          << [[error localizedDescription] UTF8String] << "\n";
                return;
            }

            available = true;
            std::cout << "[Metal] GPU backend initialized: "
                      << [[device name] UTF8String] << "\n";
        }
    }
};

MetalBackend::MetalBackend() : impl_(std::make_unique<Impl>()) {}
MetalBackend::~MetalBackend() = default;

bool MetalBackend::is_available() const { return impl_ && impl_->available; }

std::vector<bool> MetalBackend::batch_in_circumcircle(
    const std::vector<dt::Point>& vertices,
    const std::vector<std::array<int,3>>& candidates,
    const dt::Point& query_point)
{
    int n = static_cast<int>(candidates.size());
    std::vector<bool> results(n, false);

    if (n == 0) return results;

    // CPU fallback if Metal is not available
    if (!is_available()) {
        for (int i = 0; i < n; i++) {
            const auto& c = candidates[i];
            results[i] = dt::in_circumcircle(vertices[c[0]], vertices[c[1]],
                                             vertices[c[2]], query_point);
        }
        return results;
    }

    @autoreleasepool {
        auto* dev = impl_->device;

        // Prepare vertex buffer (float2 array)
        int nv = static_cast<int>(vertices.size());
        std::vector<float> verts_flat(nv * 2);
        for (int i = 0; i < nv; i++) {
            verts_flat[2*i]   = static_cast<float>(vertices[i].x);
            verts_flat[2*i+1] = static_cast<float>(vertices[i].y);
        }

        // Prepare triangle index buffer
        std::vector<int> tri_flat(n * 3);
        for (int i = 0; i < n; i++) {
            tri_flat[3*i]   = candidates[i][0];
            tri_flat[3*i+1] = candidates[i][1];
            tri_flat[3*i+2] = candidates[i][2];
        }

        // Query point
        float qp[2] = {static_cast<float>(query_point.x),
                        static_cast<float>(query_point.y)};

        // Create Metal buffers
        id<MTLBuffer> vertBuf = [dev newBufferWithBytes:verts_flat.data()
                                                length:verts_flat.size() * sizeof(float)
                                               options:MTLResourceStorageModeShared];
        id<MTLBuffer> triBuf  = [dev newBufferWithBytes:tri_flat.data()
                                                length:tri_flat.size() * sizeof(int)
                                               options:MTLResourceStorageModeShared];
        id<MTLBuffer> qpBuf   = [dev newBufferWithBytes:qp
                                                length:sizeof(qp)
                                               options:MTLResourceStorageModeShared];
        id<MTLBuffer> resBuf  = [dev newBufferWithLength:n * sizeof(int)
                                               options:MTLResourceStorageModeShared];

        // Dispatch compute
        id<MTLCommandBuffer> cmdBuf = [impl_->queue commandBuffer];
        id<MTLComputeCommandEncoder> encoder = [cmdBuf computeCommandEncoder];
        [encoder setComputePipelineState:impl_->pipeline];
        [encoder setBuffer:vertBuf offset:0 atIndex:0];
        [encoder setBuffer:triBuf  offset:0 atIndex:1];
        [encoder setBuffer:qpBuf   offset:0 atIndex:2];
        [encoder setBuffer:resBuf  offset:0 atIndex:3];

        NSUInteger threadGroupSize = std::min(
            static_cast<NSUInteger>(impl_->pipeline.maxTotalThreadsPerThreadgroup),
            static_cast<NSUInteger>(n));
        MTLSize gridSize = MTLSizeMake(n, 1, 1);
        MTLSize groupSize = MTLSizeMake(threadGroupSize, 1, 1);
        [encoder dispatchThreads:gridSize threadsPerThreadgroup:groupSize];
        [encoder endEncoding];

        [cmdBuf commit];
        [cmdBuf waitUntilCompleted];

        // Read back results
        int* resPtr = static_cast<int*>([resBuf contents]);
        for (int i = 0; i < n; i++) {
            results[i] = (resPtr[i] != 0);
        }
    }

    return results;
}
