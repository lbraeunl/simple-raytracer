#include "Renderer.hpp"

std::vector<glm::vec3> Renderer::render()
{
    std::vector<Ray> view_rays = scene.camera.generate_rays();
    std::vector<glm::vec3> pixel_values(view_rays.size());

    for (size_t i=0; i<view_rays.size();++i)
    {
        HitRecord hit = bvh.traverse_BVH(view_rays[i], false);
        if (hit.t!=INFINITY)
        {
            Ray shadow_ray(view_rays[i].at(hit.t)+hit.triangle->normal*0.0005f,scene.lights[0].position);  //rework that normal is member of HitRecord?
            HitRecord shadow_hit = bvh.traverse_BVH(shadow_ray, true);

            float intensity = 0.f;
            glm::vec3 v(0.f);
            glm::vec3 color_diffuse(0.f);
            glm::vec3 color_ambient = scene.materials[hit.triangle->mat_id].ambientColor;

            if(shadow_hit.t==INFINITY) {
                intensity = std::max(0.0f,dot(shadow_ray.direction, hit.triangle->normal));
                v = shadow_ray.direction - 2*dot(hit.triangle->normal,-shadow_ray.direction)*hit.triangle->normal;
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
            for(int c=0;c<3;++c){
                
                float spec = std::pow(std::max(0.0f, dot(-view_rays[i].direction, v)), m);

                color[c] = color_diffuse[c]*(alpha*scene.lights[0].color[c]*intensity+beta*color_ambient[c])
                +(gamma_*scene.lights[0].color[c]*intensity*std::pow(dot(-view_rays[i].direction,v),m));                   
            }

            pixel_values[i] = color;
        }
    }

    return pixel_values;
}

// Shadow rays not distance-limited

// Specular term can go NaN

// Ambient term is mathematically wrong

// Lighting normal usage inconsistent

// Back-face lighting bug