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
    Ray sample_ggx_ray(RNG& rng, const HitRecord& hit, glm::vec3 V, float roughness) const;

    glm::vec3 sample_triangle_point(RNG& rng, const Triangle& t) const;
    glm::vec3 next_event_estimation(RNG& rng, HitRecord& hit) const;

    float ggx_distribution(glm::vec3 N, glm::vec3 H, float roughness) {
        float a = roughness * roughness;
        float a2 = a * a;
        float NdotH = std::max(glm::dot(N, H), 0.0f);
        float NdotH2 = NdotH * NdotH;

        float denom = (NdotH2 * (a2 - 1.0f) + 1.0f);
        denom = std::numbers::pi_v<float> * denom * denom;

        return a2 / std::max(denom, 0.000001f);
    }

    float smith_geometry(glm::vec3 N, glm::vec3 V, glm::vec3 L, float roughness) {
        float NdotV = std::max(glm::dot(N, V), 0.0f);
        float NdotL = std::max(glm::dot(N, L), 0.0f);
        float k = (roughness * roughness) / 2.0f;
        float ggx2 = NdotV / (NdotV * (1.0f - k) + k);
        float ggx1 = NdotL / (NdotL * (1.0f - k) + k);

        return ggx1 * ggx2;
    }

    glm::vec3 schlick_fresnel(float cos_theta, glm::vec3 F0) {
        return F0 + (1.0f - F0) * std::pow(std::clamp(1.0f - cos_theta, 0.0f, 1.0f), 5.0f);
    }
};