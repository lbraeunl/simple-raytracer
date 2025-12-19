#include "BVH.hpp"
#include <stack>



int BVH::build_BVH(int start, int end)
{
    bool isLeaf = (end - start) <= leafsize;
    nodes.emplace_back(start, end, isLeaf);
    int nodeIndex = nodes.size() - 1;

    nodes[nodeIndex].box = compute_box(start, end);

    if (isLeaf)
    {
        return nodeIndex;
    }

    int axis = nodes[nodeIndex].box.longest_axis();
    size_t mid = (start+end)/2;

    std::nth_element(indices.begin()+start, indices.begin()+mid, indices.begin()+end, [axis, this](int a, int b) {
        return (*triangles)[a].centroid[axis] < (*triangles)[b].centroid[axis];
    });

    nodes[nodeIndex].left = build_BVH(start, mid);
    nodes[nodeIndex].right = build_BVH(mid, end);

    return nodeIndex;
}


HitRecord BVH::traverse_BVH(const Ray& ray, bool any_hit) const
{
    HitRecord best_hit(INFINITY, nullptr);
    std::stack<int> stack;
    stack.push(root);

    //LOG(heatMap.intersects = 0;);

    while (!stack.empty()) {
        int nodeIndex = stack.top();
        stack.pop();
        const BVHNode& node = nodes[nodeIndex];

        //LOG(heatMap.intersects += 1;);

        if (!(node.box.box_intersect(ray)))
            continue;

        if (node.isLeaf) {
            for (int i = 0; i < (node.end-node.start); i++)
            {
                const Triangle& t = (*triangles)[indices[node.start+i]];
                HitRecord current_hit;
                if (t.triangle_intersect(ray, current_hit) && current_hit.t < best_hit.t) {
                    best_hit = current_hit;
                    if (any_hit) 
                        return best_hit;
                }
            }

        } else {
            stack.push(node.left);
            stack.push(node.right);
        }
    }
    return best_hit;
}


AABB BVH::compute_box(int start, int end)
{
    AABB box;
    for (int i = start; i < end; ++i)
        box.expand((*triangles)[indices[i]]);
    return box;
}

