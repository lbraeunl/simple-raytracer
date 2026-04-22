#include "Renderer.hpp"
#include "logger.hpp"
#include <thread>
#include <atomic>

std::vector<glm::vec3> Renderer::single_render() 
{ 
    const int width  = scene.camera.resolution.x;
    const int height = scene.camera.resolution.y;
    std::vector<glm::vec3> framebuffer(width * height);
    LOG(heat_map.heatpixels.resize(width * height));

    for (size_t py = 0; py < height; ++py) {
        for (size_t px = 0; px < width; ++px) {

            glm::vec3 color = trace(px, py);

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

                    glm::vec3 color = trace(px, py);

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
        return scene.materials[hit.triangle->mat_id].diffuseColor;
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
    int hit_count = 0;

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