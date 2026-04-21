#pragma once
#include <numeric>
#include "Geometry.hpp"



struct Bucket {
    // int start = -1;
    // int end = -1;
    int count = 0;
    AABB box{};
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

        int end = (int)tris.size();
        AABB root_box = compute_box(0, end);
        root = build_BVH(0, end, root_box);
    }

    HitRecord traverse_BVH(const Ray& ray, float t_max) const;

private:
    int build_BVH(int start, int end, AABB parent_box, int nBuckets = 12);
    int build_simple_BVH(int start, int end);
    AABB compute_box(int start, int end);
    AABB compute_box_fast(int start, int end, int axis);   
};