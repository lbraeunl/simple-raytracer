#pragma once
#include "Scene.hpp"
#include "BVH.hpp"
#include "utils.hpp"


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

    const int sample_count = 10;
    const int bounces = 1;
    const float diffuse_probability = 0.5f;

    Renderer(const Scene& scene_, const BVH& bvh_) : scene(scene_), bvh(bvh_) {}

    std::vector<glm::vec3> render();
    std::vector<glm::vec3> single_render();

private:
    glm::vec3 trace(int px, int py) const;
    glm::vec3 stochastic_trace(int px, int py) const;
    glm::vec3 get_material_color(const HitRecord& hit) const;
    glm::vec3 compute_radiance(Ray ray, HitRecord hit, RNG& rng) const;

    glm::vec3 sample_triangle_point(RNG& rng, const Triangle& t) const;
    glm::vec3 next_event_estimation(RNG& rng, HitRecord& hit, const glm::vec3& wo, const glm::vec3& base_color) const;
};