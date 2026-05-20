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


glm::vec3 Renderer::get_material_color(HitRecord& hit) const
{
    int tex_id = scene.materials[hit.triangle->mat_id].diffuseTex;
    if(tex_id >=0)
    {   
        float w = 1.0f - hit.u - hit.v;
        glm::vec2 uv_values = w * hit.triangle->uv[0] + hit.u * hit.triangle->uv[1] + hit.v * hit.triangle->uv[2];
        const Texture& tex = scene.textures[tex_id];

        float u = uv_values.x - std::floor(uv_values.x);
        float v = uv_values.y - std::floor(uv_values.y);

        int x = int(u * (tex.width - 1));
        int y = int((v) * (tex.height - 1));

        int tex_idx = (y * tex.width + x) * tex.channels;

        return {tex.data[tex_idx + 0],tex.data[tex_idx + 1],tex.data[tex_idx + 2]};
    }
    else
    {
        return scene.materials[hit.triangle->mat_id].albedo;
    }
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


glm::vec3 Renderer::compute_radiance(Ray& ray, HitRecord& hit, RNG& rng, int depth) const
{
    glm::vec3 radiance = scene.materials[hit.triangle->mat_id].emissive;
    if (radiance.x != 0.0f || radiance.y != 0.0f || radiance.z != 0.0f) return radiance;
    
    glm::vec3 rho = scene.materials[hit.triangle->mat_id].albedo;
    float e = scene.materials[hit.triangle->mat_id].roughness;
    glm::vec3 n = hit.normal;

    for (int j=0; j < sample_count; j++)
    {
        radiance += (1.0f/sample_count)*next_event_estimation(rng, hit);

        //Ray secondary_ray = sample_cosine_weighted_ray(rng, hit);
        Ray secondary_ray = sample_glossy_ray(rng, hit,ray.direction,e);

        HitRecord secondary_hit = bvh.traverse_BVH(secondary_ray, INFINITY);

        if (secondary_hit.t == INFINITY)
            continue;
        
        glm::vec3 l_i = scene.materials[secondary_hit.triangle->mat_id].emissive;

        if (depth < bounces)
        {
            l_i = compute_radiance(secondary_ray, secondary_hit, rng, depth+1);
        }

        glm::vec3 wi = secondary_ray.direction;
        glm::vec3 r =glm::reflect(-ray.direction, n);
        float cos_theta = glm::max(glm::dot(n, wi), 0.0f);
        float phong_cos = glm::max(glm::dot(r, wi), 0.0f);

        glm::vec3 brdf = rho *((e + 2.0f) /(2.0f * std::numbers::pi_v<float>)) *std::pow(phong_cos, e);

        float pdf = ((e + 1.0f) /(2.0f * std::numbers::pi_v<float>)) *std::pow(phong_cos, e);

        radiance += (1.0f / sample_count)*brdf*l_i*cos_theta /pdf;

        //radiance += (1.0f/sample_count)*rho*l_i;
    }

    return radiance;
}


Ray Renderer::sample_cosine_weighted_ray(RNG& rng, HitRecord& hit) const
{
    float xi1 = rng.next_float();
    float xi2 = rng.next_float();
    float phi = 2.0*std::numbers::pi*xi2;
    float r = sqrt(xi1);

    glm::vec3 local_p2(r*cos(phi),r*sin(phi),sqrt(1.0-xi1));
    glm::quat q = glm::rotation(glm::vec3(0.f,0.f,1.0f),hit.normal);

    Ray ray(hit.position + hit.normal*0.00001f, hit.position + q*local_p2);

    return ray;
}

Ray Renderer::sample_glossy_ray(RNG& rng, const HitRecord& hit, const glm::vec3& wo, float e) const
{
    glm::vec3 r = glm::reflect(-wo, hit.normal);

    float xi1 = rng.next_float();
    float xi2 = rng.next_float();
    float cos_theta = std::pow(xi1, 1.0f / (e + 1.0f));
    float sin_theta = std::sqrt(std::max(0.0f, 1.0f - cos_theta * cos_theta));
    float phi = 2.0f*std::numbers::pi_v<float>*xi2;

    glm::vec3 local_p2(sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta);
    glm::quat q = glm::rotation(glm::vec3(0.f, 0.f, 1.f), r);

    Ray ray(hit.position + r * 0.00001f, hit.position + q*local_p2);

    return ray;
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

glm::vec3 Renderer::next_event_estimation(RNG& rng, HitRecord& hit) const
{
    glm::vec3 l_i(0.f);
    if (scene.lighted_triangles.empty()) return l_i;

    int light_count = scene.lighted_triangles.size();
    int n = int(rng.next_float()*light_count);
    size_t idx = scene.lighted_triangles[n];
    const Triangle& t = scene.triangles[idx];

    glm::vec3 light_pos = sample_triangle_point(rng, t);

    Ray shadow_ray(hit.position + hit.normal*0.0001f, light_pos);
    float light_distance = glm::distance(light_pos, hit.position+ hit.normal*0.001f);

    HitRecord shadow_hit = bvh.traverse_BVH(shadow_ray, light_distance);

    if(shadow_hit.t < light_distance) return l_i;
    
    glm::vec3 l_e = scene.materials[t.mat_id].emissive;
    glm::vec3 rho = scene.materials[hit.triangle->mat_id].albedo;
    glm::vec3 w_i = glm::normalize(light_pos - hit.position);
    float cos_surf = std::max(glm::dot(hit.normal, w_i), 0.0f);
    float cos_light = std::max(glm::dot(t.normal, -w_i), 0.0f);

    l_i = l_e*rho/std::numbers::pi_v<float>*(cos_surf*cos_light)/(light_distance*light_distance)*(light_count * t.area);
    return l_i;
}

glm::vec3 Renderer::stochastic_trace(int px, int py) const
{
    Ray ray = scene.camera.generate_ray(px, py);

    HitRecord hit = bvh.traverse_BVH(ray, INFINITY);

    if (hit.t == INFINITY)
        return glm::vec3(0.0f);

    thread_local RNG rng(((px * 1597334677U) ^ (py * 3812015801U))+1);

    glm::vec3 radiance = compute_radiance(ray, hit, rng, 1);

    return radiance;
}