#include "Renderer.hpp"

std::vector<glm::vec3> Renderer::render()
{
    std::vector<Ray> view_rays = scene.camera.generate_rays();
    std::vector<glm::vec3> pixel_values(view_rays.size());

    for (size_t i=0; i<view_rays.size();++i)
    {
        HitRecord hit = bvh.traverse_BVH(view_rays[i], false);

        if (hit.t == INFINITY)
            continue;

        glm::vec3 normal = hit.triangle->interpolate_normal(hit.u, hit.v);

        Ray shadow_ray(view_rays[i].at(hit.t)+normal*0.0005f,scene.lights[0].position);
        HitRecord shadow_hit = bvh.traverse_BVH(shadow_ray, true);

        float diff_intensity = 0.f;
        float spec_intensity = 0.f;

        glm::vec3 color_diffuse(0.f);
        glm::vec3 color_ambient = scene.materials[hit.triangle->mat_id].ambientColor;

        if(shadow_hit.t==INFINITY) { //change this later, so it only counts until the light source.
            glm::vec3 r = -shadow_ray.direction + 2*dot(normal,shadow_ray.direction)*normal;
            diff_intensity = std::max(0.0f,dot(shadow_ray.direction, normal));
            spec_intensity = std::pow(std::max(0.0f, dot(-view_rays[i].direction, r)), m);
        }

        int tex_id = scene.materials[hit.triangle->mat_id].diffuseTex;
        if(tex_id >=0)
        {
            glm::vec2 uv_values = hit.u*hit.triangle->uv[0] + hit.v*hit.triangle->uv[1] + (1-hit.u-hit.v) * hit.triangle->uv[2];
            Texture tex = scene.textures[tex_id];

            int x = glm::clamp(int(uv_values.x * (tex.width  - 1)), 0, tex.width  - 1);
            int y = glm::clamp(int(uv_values.y * (tex.height - 1)), 0, tex.height - 1);
            int tex_idx = (y * tex.width + x) * tex.channels;

            color_diffuse = {tex.data[tex_idx + 0],tex.data[tex_idx + 1],tex.data[tex_idx + 2]};
        }
        else
        {
            color_diffuse = scene.materials[hit.triangle->mat_id].diffuseColor;
        }

        glm::vec3 color;
        for(int c=0;c<3;++c)
        {
            color[c] =
            alpha*color_diffuse[c] +
            beta*scene.lights[0].color[c]*diff_intensity*color_diffuse[c] +
            gamma_*scene.lights[0].color[c]*spec_intensity;
        }
        pixel_values[i] = color;
    }

    return pixel_values;
}