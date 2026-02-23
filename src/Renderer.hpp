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

    float alpha = 0.3;
    float beta = 5;
    float gamma_ = 5;
    float m = 10;

    Renderer(const Scene& scene_, const BVH& bvh_) : scene(scene_), bvh(bvh_) {}

    std::vector<glm::vec3> render();
    std::vector<glm::vec3> single_render();

private:
    HitRecord shadow_ray();
    glm::vec3 trace(Ray ray) const;
    glm::vec3 get_material_color(HitRecord& hit) const;
};