#include "Material.hpp"
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include <numbers>
#include <algorithm>

glm::vec3 Material::eval(const glm::vec3& wo, const glm::vec3& wi, const glm::vec3& normal) const
{
    return albedo / std::numbers::pi_v<float>;
}

glm::vec3 Material::diffuse_sample(RNG& rng, const glm::vec3& wo, const glm::vec3& normal, glm::vec3& wi, float& pdf) const
{
    wi = BRDF::sample_cosine_weighted(rng, normal);

    float cos_theta = glm::max(glm::dot(normal, wi), 0.0f);
    pdf = cos_theta / std::numbers::pi_v<float>;

    return albedo / std::numbers::pi_v<float>;
}

glm::vec3 Material::sample(RNG& rng, const glm::vec3& wo, const glm::vec3& normal, glm::vec3& wi, float& pdf) const
{
    glm::vec3 r = -glm::reflect(wo, normal);
    
    wi = BRDF::sample_glossy(rng, r, roughness);

    float cos_theta = glm::max(glm::dot(normal, wi), 0.0f);
    float phong_cos = glm::max(glm::dot(r, wi), 0.0f);

    pdf = ((roughness + 1.0f) / (2.0f * std::numbers::pi_v<float>)) * std::pow(phong_cos, roughness);

    glm::vec3 brdf = albedo * ((roughness + 2.0f) / (2.0f * std::numbers::pi_v<float>)) * std::pow(phong_cos, roughness);

    return brdf;
}

namespace BRDF {

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

glm::vec3 sample_cosine_weighted(RNG& rng, const glm::vec3& normal) {
    float xi1 = rng.next_float();
    float xi2 = rng.next_float();
    float phi = 2.0f * std::numbers::pi_v<float> * xi2;
    float r = std::sqrt(xi1);

    glm::vec3 local_p2(r * std::cos(phi), r * std::sin(phi), std::sqrt(1.0f - xi1));
    glm::quat q = glm::rotation(glm::vec3(0.f, 0.f, 1.0f), normal);

    return glm::normalize(q * local_p2);
}

glm::vec3 sample_glossy(RNG& rng, const glm::vec3& r, float e) {
    float xi1 = rng.next_float();
    float xi2 = rng.next_float();
    float cos_theta = std::pow(xi1, 1.0f / (e + 1.0f));
    float sin_theta = std::sqrt(std::max(0.0f, 1.0f - cos_theta * cos_theta));
    float phi = 2.0f * std::numbers::pi_v<float> * xi2;

    glm::vec3 local_p2(sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta);
    glm::quat q = glm::rotation(glm::vec3(0.f, 0.f, 1.f), r);

    return glm::normalize(q * local_p2);
}

glm::vec3 sample_ggx(RNG& rng, const glm::vec3& normal, glm::vec3 V, float roughness) {
    float a = roughness * roughness;
    float xi1 = rng.next_float();
    float xi2 = rng.next_float();

    float phi = 2.0f * std::numbers::pi_v<float> * xi1;
    float cos_theta = std::sqrt((1.0f - xi2) / (1.0f + (a * a - 1.0f) * xi2));
    float sin_theta = std::sqrt(1.0f - cos_theta * cos_theta);

    glm::vec3 H_tangent(sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta);

    glm::quat q = glm::rotation(glm::vec3(0.f, 0.f, 1.f), normal);
    glm::vec3 H = glm::normalize(q * H_tangent);

    return glm::normalize(glm::reflect(-V, H));
}

} // namespace BRDF
