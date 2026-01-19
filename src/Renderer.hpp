#pragma once
#include "Scene.hpp"
#include "BVH.hpp"


class Renderer {
public:
    const Scene& scene;
    const BVH& bvh;

    float alpha = 1;
    float beta = 6;
    float gamma_ = 5;
    float m = 10;

    Renderer(const Scene& scene_, const BVH& bvh_) : scene(scene_), bvh(bvh_) {}

    std::vector<glm::vec3> render();
    glm::vec3 illumination(HitRecord hit, HitRecord shadow_hit);

private:
    HitRecord shadow_ray();
};