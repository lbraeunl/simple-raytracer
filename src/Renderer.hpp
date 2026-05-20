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

    int sample_count = 40;
    int bounces = 1;

    Renderer(const Scene& scene_, const BVH& bvh_) : scene(scene_), bvh(bvh_) {}

    std::vector<glm::vec3> render();
    std::vector<glm::vec3> single_render();

private:
    glm::vec3 trace(int px, int py) const;
    glm::vec3 stochastic_trace(int px, int py) const;
    glm::vec3 get_material_color(HitRecord& hit) const;
    glm::vec3 compute_radiance(Ray& ray, HitRecord& hit, RNG& rng, int depth) const;

    Ray sample_cosine_weighted_ray(RNG& rng, HitRecord& hit) const;
    Ray sample_glossy_ray(RNG& rng, const HitRecord& hit, const glm::vec3& wo, float e) const;

    glm::vec3 sample_triangle_point(RNG& rng, const Triangle& t) const;
    glm::vec3 next_event_estimation(RNG& rng, HitRecord& hit) const;
};