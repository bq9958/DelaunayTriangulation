#pragma once

#include <unordered_map>
#include <cstdint>
#include <utility>

// Replaces the custom HashTable from the C version.
// Maps directed edge (v1,v2) -> (triangle_id, local_edge_index).
class EdgeMap {
public:
    struct Info {
        int tri;
        int edge;
    };

    void add(int v1, int v2, int tri, int edge) {
        map_[key(v1, v2)] = {tri, edge};
    }

    // Directed lookup: exact (v1,v2) match.
    Info* find_directed(int v1, int v2) {
        auto it = map_.find(key(v1, v2));
        return it != map_.end() ? &it->second : nullptr;
    }

    // Undirected lookup: tries (v1,v2) then (v2,v1).
    Info* find(int v1, int v2) {
        if (auto* p = find_directed(v1, v2)) return p;
        return find_directed(v2, v1);
    }

    void remove(int v1, int v2) {
        map_.erase(key(v1, v2));
    }

    void clear() { map_.clear(); }
    size_t size() const { return map_.size(); }

private:
    static uint64_t key(int v1, int v2) {
        return (static_cast<uint64_t>(static_cast<uint32_t>(v1)) << 32)
             |  static_cast<uint64_t>(static_cast<uint32_t>(v2));
    }

    std::unordered_map<uint64_t, Info> map_;
};
