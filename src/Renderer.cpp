#include "Renderer.hpp"
#include "logger.hpp"
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include <thread>
#include <atomic>
#include <numbers>

std::vector<glm::vec3> Renderer::single_render() 
{ 
    const int width  = scene.camera.resolution.x;
    const int height = scene.camera.resolution.y;
    std::vector<glm::vec3> framebuffer(width * height);
    LOG(heat_map.heatpixels.resize(width * height));

    for (size_t py = 0; py < height; ++py) {
        for (size_t px = 0; px < width; ++px) {

            glm::vec3 color = stochastic_trace(px, py);

            framebuffer[py * width + px] = color;
        }
    }
    
    return framebuffer; 
}

    
std::vector<glm::vec3> Renderer::render()
{
    const int width  = scene.camera.resolution.x;
    const int height = scene.camera.resolution.y;
    const int tile_size = 32;

    std::vector<glm::vec3> framebuffer(width * height);

    std::vector<Tile> tiles;
    for (int y = 0; y < height; y += tile_size) 
    {
        for (int x = 0; x < width; x += tile_size) 
        {
            tiles.push_back({x,y,std::min(x + tile_size, width),std::min(y + tile_size, height)});
        }
    }

    std::atomic<size_t> nextTile{0};
    const unsigned thread_count = std::thread::hardware_concurrency();

    auto worker = [&]() {
        while (true) {
            size_t tileIndex = nextTile.fetch_add(1);
            if (tileIndex >= tiles.size())
                break;

            const Tile& tile = tiles[tileIndex];

            for (size_t py = tile.y0; py < tile.y1; ++py) {
                for (size_t px = tile.x0; px < tile.x1; ++px) {

                    glm::vec3 color = stochastic_trace(px, py);

                    framebuffer[py * width + px] = color;
                }
            }
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(thread_count);

    for (unsigned i = 0; i < thread_count; ++i)
        threads.emplace_back(worker);

    for (auto& t : threads)
        t.join();

    return framebuffer;
}


glm::vec3 Renderer::get_material_color(const HitRecord& hit) const
{
    const Material& material = scene.materials[hit.triangle->mat_id];
    if (material.diffuseTex < 0)
        return material.albedo;

    float w = 1.0f - hit.u - hit.v;
    glm::vec2 uv = w * hit.triangle->uv[0] + hit.u * hit.triangle->uv[1] + hit.v * hit.triangle->uv[2];

    return material.albedo * scene.textures[material.diffuseTex].sample(uv);
}


glm::vec3 Renderer::trace(int px, int py) const
{
    Ray ray = scene.camera.generate_ray(px, py);

    LOG(INT_COUNT += heat_map.intersects);
    LOG(heat_map.intersects = 0);
    HitRecord hit = bvh.traverse_BVH(ray, INFINITY);
    LOG(heat_map.heatpixels[py * scene.camera.resolution.x + px] = heat_map.intersects);

    if (hit.t == INFINITY)
        return glm::vec3(0.0f);

    glm::vec3 mat_color = get_material_color(hit);
    glm::vec3 normal = hit.triangle->interpolate_normal(hit.u, hit.v);

    glm::vec3 total_diffuse(0.f);
    glm::vec3 total_specular(0.f);

    for (const LightSource& l : scene.lights)
    {
        auto light_points = l.get_lightpoints();
        int sample_count = light_points.size();

        for (const glm::vec3& lp : light_points)
        {
            glm::vec3 hit_position = ray.at(hit.t);            
            Ray shadow_ray(hit_position+normal*0.0005f,lp);
            float light_distance = glm::distance(lp, hit_position);

            HitRecord shadow_hit = bvh.traverse_BVH(shadow_ray, light_distance);

            float diffuse_factor = 0.f;
            float specular_factor = 0.f;

            if(shadow_hit.t >= light_distance) {
                glm::vec3 r = -shadow_ray.direction + 2*dot(normal,shadow_ray.direction)*normal;
                diffuse_factor = std::max(0.0f,dot(shadow_ray.direction, normal));
                specular_factor = std::pow(std::max(0.0f, dot(-ray.direction, r)), m);
            }

            total_diffuse += l.color*mat_color*diffuse_factor;
            total_specular += l.color*specular_factor;
        }

        total_diffuse /= sample_count;
        total_specular /= sample_count;
    }

    glm::vec3 color = alpha*mat_color + beta*total_diffuse + gamma_*total_specular;
    
    return color;
}


glm::vec3 Renderer::compute_radiance(Ray ray, HitRecord hit, RNG& rng) const
{
    glm::vec3 radiance(0.0f);
    glm::vec3 throughput(1.0f);
    float prev_pdf_brdf = 1.0f;
    bool is_specular_bounce = true;

    for (int depth = 0; depth <= bounces; depth++)
    {
        glm::vec3 emission = scene.materials[hit.triangle->mat_id].emissive;
        if (emission.x != 0.0f || emission.y != 0.0f || emission.z != 0.0f) {
            if (!hit.back_face) {
                if (is_specular_bounce) {
                    radiance += throughput * emission;
                } else {
                    float light_distance = hit.t;
                    float cos_light = std::max(glm::dot(hit.triangle->normal, -ray.direction), 0.0f);
                    int light_count = scene.lighted_triangles.size();
                    float light_area = hit.triangle->area;

                    float pdf_light = 0.0f;
                    if (light_area > 0.0f && light_count > 0) {
                        pdf_light = (light_distance * light_distance) / (cos_light * light_area * light_count);
                    }

                    float w_brdf = 1.0f;
                    if (pdf_light > 0.0f) {
                        w_brdf = (prev_pdf_brdf * prev_pdf_brdf) / (prev_pdf_brdf * prev_pdf_brdf + pdf_light * pdf_light);
                    }
                    radiance += throughput * emission * w_brdf;
                }
            }
            break;
        }

        if (depth == bounces) break;

        glm::vec3 wo = -ray.direction;
        glm::vec3 base_color = get_material_color(hit);
        radiance += throughput * next_event_estimation(rng, hit, wo, base_color);

        glm::vec3 wi;
        float pdf_brdf;
        const Material& material = scene.materials[hit.triangle->mat_id];
        glm::vec3 brdf = material.sample(rng, wo, hit.normal, base_color, wi, pdf_brdf);

        if (pdf_brdf <= 0.0f) break;

        float cos_theta = glm::max(glm::dot(hit.normal, wi), 0.0f);
        
        throughput *= brdf * cos_theta / pdf_brdf;
        prev_pdf_brdf = pdf_brdf;
        is_specular_bounce = false;

        float max_throughput = std::max({throughput.x, throughput.y, throughput.z});
        if (depth > 2) {
            float survival_prob = std::min(1.0f, max_throughput);
            if (rng.next_float() > survival_prob) {
                break;
            }
            throughput /= survival_prob;
        }

        ray = Ray(hit.position + hit.normal * 0.0001f, hit.position + wi);
        hit = bvh.traverse_BVH(ray, INFINITY);

        if (hit.t == INFINITY) break;
    }

    return radiance;
}


glm::vec3 Renderer::sample_triangle_point(RNG& rng, const Triangle& t) const
{
    float xi1 = rng.next_float();
    float xi2 = rng.next_float();
    float su = sqrt(1.0f - xi1);

    float u = 1.0f - su;
    float v = xi2 * su;
    float w = 1.0f - u - v;

    return u * t.v[0] + v * t.v[1] + w * t.v[2];
}


glm::vec3 Renderer::next_event_estimation(RNG& rng, HitRecord& hit, const glm::vec3& wo, const glm::vec3& base_color) const
{
    glm::vec3 radiance(0.f);
    if (scene.lighted_triangles.empty()) return radiance;

    int light_count = scene.lighted_triangles.size();
    int n = int(rng.next_float() * light_count);
    size_t idx = scene.lighted_triangles[n];
    const Triangle& t = scene.triangles[idx];

    glm::vec3 light_pos = sample_triangle_point(rng, t);

    Ray shadow_ray(hit.position + hit.normal * 0.0001f, light_pos);
    float light_distance = glm::distance(light_pos, hit.position + hit.normal * 0.001f);

    HitRecord shadow_hit = bvh.traverse_BVH(shadow_ray, light_distance);

    if (shadow_hit.t < light_distance) return radiance;

    glm::vec3 l_e = scene.materials[t.mat_id].emissive;
    glm::vec3 w_i = glm::normalize(light_pos - hit.position);
    float cos_surf = std::max(glm::dot(hit.normal, w_i), 0.0f);
    float cos_light = std::max(glm::dot(t.normal, -w_i), 0.0f);

    if (cos_light <= 0.0f) return radiance;

    float pdf_light = (light_distance * light_distance) / (cos_light * t.area * light_count);
    const Material& material = scene.materials[hit.triangle->mat_id];
    float pdf_brdf = material.compute_pdf(wo, w_i, hit.normal);
    float w_light = (pdf_light * pdf_light) / (pdf_light * pdf_light + pdf_brdf * pdf_brdf);

    glm::vec3 brdf = material.eval(wo, w_i, hit.normal, base_color);
    radiance = l_e * brdf * cos_surf / pdf_light * w_light;

    return radiance;
}


glm::vec3 Renderer::stochastic_trace(int px, int py) const
{
    glm::vec3 total_radiance(0.0f);
    thread_local RNG rng(((px * 1597334677U) ^ (py * 3812015801U))+1);

    for (int i = 0; i < sample_count; i++) {
        float jx = px;//+ rng.next_float();
        float jy = py; //+ rng.next_float();
        Ray ray = scene.camera.generate_ray(jx, jy);

        HitRecord hit = bvh.traverse_BVH(ray, INFINITY);

        if (hit.t != INFINITY) {
            total_radiance += compute_radiance(ray, hit, rng);
        }
    }

    return total_radiance / float(sample_count);
}