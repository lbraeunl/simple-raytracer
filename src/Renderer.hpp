#pragma once
#include "Scene.hpp"
#include "BVH.hpp"

struct Tile 
{ 
    int x0, y0, x1, y1;
};


class Renderer {
public:
    const Scene& scene;
    const BVH& bvh;

    float alpha = 0.6;
    float beta = 3;
    float gamma_ = 3;
    float m = 10;

    Renderer(const Scene& scene_, const BVH& bvh_) : scene(scene_), bvh(bvh_) {}

    std::vector<glm::vec3> render();
    std::vector<glm::vec3> single_render();

private:
    glm::vec3 trace(int px, int py) const;
    glm::vec3 get_material_color(HitRecord& hit) const;
};