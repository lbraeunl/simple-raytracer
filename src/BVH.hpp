#pragma once
#include <numeric>
#include "Geometry.hpp"



struct Bucket {
    int count = 0;
    AABB box;
};


struct BVHNode {
    AABB box;
    int start;
    int end;
    int left = -1;
    int right = -1;
    bool isLeaf;

    BVHNode(int start, int end, bool isLeaf=false): start(start), end(end), isLeaf(isLeaf) {}
};


class BVH {
public:
    const std::vector<Triangle>* triangles;
    std::vector<int> indices;
    std::vector<BVHNode> nodes;
    int root;
    int leafsize = 4;

    BVH(const std::vector<Triangle>& tris)
    {
        triangles = &tris;
        indices.resize(tris.size());
        std::iota(indices.begin(), indices.end(), 0);
        nodes.reserve(tris.size() * 2);
        root = build_BVH(0, (int)tris.size());
    }

    HitRecord traverse_BVH(const Ray& ray, bool any_hit);

private:
    int build_BVH(int start, int end);
    AABB compute_box(int start, int end);
};