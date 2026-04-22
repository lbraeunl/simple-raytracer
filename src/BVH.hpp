#pragma once
#include <numeric>
#include "Geometry.hpp"



struct Bucket {
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
    std::vector<AABB> boxes;
    int root;
    int leafsize = 4;

    BVH(const std::vector<Triangle>& tris)
    {
        triangles = &tris;
        int n = (int)tris.size();

        indices.resize(n);
        std::iota(indices.begin(), indices.end(), 0);
        nodes.reserve(n*2);

        AABB root_box;
        for (int i = 0; i < n; i++) {
            boxes.emplace_back(tris[i]);
            root_box.expand(boxes[i]);
        }       

        root = build_simple_BVH(0, n, root_box);

        boxes.clear();
    }

    HitRecord traverse_BVH(const Ray& ray, float t_max) const;

private:
    int build_BVH(int start, int end, AABB parent_box, int nBuckets = 12);
    int build_simple_BVH(int start, int end, AABB parent_box);
    int build_random_BVH(int start, int end, AABB parent_box);

    AABB compute_box(int start, int end);
};