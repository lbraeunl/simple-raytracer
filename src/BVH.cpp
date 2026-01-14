#include "BVH.hpp"
#include <stack>


int BVH::build_BVH(int start, int end, int nBuckets)
{   
    bool isLeaf = (end - start) <= leafsize;
    nodes.emplace_back(start, end, isLeaf);
    int nodeIndex = nodes.size() - 1;

    nodes[nodeIndex].box = compute_box(start, end);

    if (isLeaf) return nodeIndex;

    int axis = nodes[nodeIndex].box.longest_axis();

    std::sort(indices.begin()+start, indices.begin()+end, [axis, this](int a, int b) {
        return (*triangles)[a].centroid[axis] < (*triangles)[b].centroid[axis];
    });

    int best_split = 0;
    float cmin = (*triangles)[indices[start]].centroid[axis];
    float cmax = (*triangles)[indices[end-1]].centroid[axis];
    float extent = cmax - cmin;

    if (extent==0){
        best_split = (start+end)/2;
    }
    else
    {
        std::vector<Bucket> buckets(nBuckets);
        for(size_t i = start; i<end; ++i)
        {
            Triangle t = (*triangles)[indices[i]];
            int b = nBuckets*(t.centroid[axis]-cmin) / extent;
            if (b == nBuckets) b = nBuckets - 1;
            buckets[b].count++;
            buckets[b].box.expand(t);
        }

        float minCost = INFINITY;
        int count = 0;
        for(int splits = 1;splits<buckets.size();++splits)
        {
            AABB left_box = buckets[0].box;
            for (int i = 1;i<splits;++i)
            {
                left_box = AABB(left_box,buckets[i].box);
            }

            AABB right_box = buckets[splits].box;
            for (int i = splits+1;i<buckets.size();++i)
            {
                right_box = AABB(right_box,buckets[i].box);
            }

            count += buckets[splits-1].count;
            float cost = count*left_box.surface_area()+((end-start)-count)*right_box.surface_area();
            if(cost < minCost){
                minCost = cost;
                best_split = start + count;
            }          
        }
    }

    nodes[nodeIndex].left = build_BVH(start, best_split, nBuckets);
    nodes[nodeIndex].right = build_BVH(best_split, end, nBuckets);

    return nodeIndex;
}


int BVH::build_simple_BVH(int start, int end)
{
    bool isLeaf = (end - start) <= leafsize;
    nodes.emplace_back(start, end, isLeaf);
    int nodeIndex = nodes.size() - 1;

    nodes[nodeIndex].box = compute_box(start, end);

    if (isLeaf) return nodeIndex;

    int axis = nodes[nodeIndex].box.longest_axis();
    size_t mid = (start+end)/2;

    std::nth_element(indices.begin()+start, indices.begin()+mid, indices.begin()+end, [axis, this](int a, int b) {
        return (*triangles)[a].centroid[axis] < (*triangles)[b].centroid[axis];
    });

    nodes[nodeIndex].left = build_simple_BVH(start, mid);
    nodes[nodeIndex].right = build_simple_BVH(mid, end);

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
            for (int i = node.start; i < node.end; ++i)
            {
                const Triangle& t = (*triangles)[indices[i]];
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