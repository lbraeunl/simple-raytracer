#pragma once
#include <glm/glm.hpp>
#include "utils.hpp"

struct Material 
{
    glm::vec3 albedo{1.0f};
    glm::vec3 emissive{0.0f};
    float metallic = 0.0f;
    float roughness = 0.5f;
    float ior = 1.0f;

    int ambientTex = -1;
    int diffuseTex = -1;
    int specularTex = -1;

    glm::vec3 eval(const glm::vec3& wo, const glm::vec3& wi, const glm::vec3& normal) const;
    float compute_pdf(const glm::vec3& wo, const glm::vec3& wi, const glm::vec3& normal, float p_diffuse) const;
    glm::vec3 sample(RNG& rng, const glm::vec3& wo, const glm::vec3& normal, glm::vec3& wi, float& pdf) const;
    glm::vec3 diffuse_sample(RNG& rng, const glm::vec3& wo, const glm::vec3& normal, glm::vec3& wi, float& pdf) const;
};

namespace BRDF {
    float ggx_distribution(glm::vec3 N, glm::vec3 H, float roughness);
    float smith_geometry(glm::vec3 N, glm::vec3 V, glm::vec3 L, float roughness);
    glm::vec3 schlick_fresnel(float cos_theta, glm::vec3 F0);

    glm::vec3 sample_cosine_weighted(RNG& rng, const glm::vec3& normal);
    glm::vec3 sample_glossy(RNG& rng, const glm::vec3& r, float e);
    glm::vec3 sample_ggx(RNG& rng, const glm::vec3& normal, glm::vec3 V, float roughness);
}
