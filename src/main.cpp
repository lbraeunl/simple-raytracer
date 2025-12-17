#include <SFML/Graphics.hpp>
#include "tiny_obj_loader.h"
#include <iostream>
#include <glm/vec3.hpp>
#include "Geometry.hpp"
#define GLM_ENABLE_EXPERIMENTAL
#define STB_IMAGE_IMPLEMENTATION
#include <glm/gtx/norm.hpp>
#include <algorithm>
#include <stack>
#include <chrono>
#include <tinycolormap.hpp>
#include <optional>
#include "Model.hpp"
#include "Scene.hpp"
#include "BVH.hpp"

using namespace glm;

constexpr bool HEATMAP = false;
#define LOG(code) \
    do { if constexpr (HEATMAP) { code; } } while(0)
struct HeatMap {
    uint8_t intersects;
    std::vector<uint8_t> heatpixels = {};
};
inline thread_local HeatMap heatMap;


std::vector<uint8_t> generateHeatmap(const std::vector<uint8_t>& pixels) {
    std::vector<uint8_t> heatmap;
    heatmap.reserve(pixels.size()*4);

    uint8_t minPixel = *std::min_element(pixels.begin(), pixels.end());
    uint8_t maxPixel = *std::max_element(pixels.begin(), pixels.end());
    float range = (maxPixel > minPixel) ? float(maxPixel - minPixel) : 1.0f;

    for (uint8_t val : pixels) {
        float normalized = (float(val) - float(minPixel)) / range;
        tinycolormap::Color c = tinycolormap::GetColor(normalized, tinycolormap::ColormapType::Turbo);

        heatmap.push_back(static_cast<uint8_t>(c.r() * 255.0f));
        heatmap.push_back(static_cast<uint8_t>(c.g() * 255.0f));
        heatmap.push_back(static_cast<uint8_t>(c.b() * 255.0f));
        heatmap.push_back(255);
    }
    return heatmap;
}

// BVHNode* buildBVH_CM(const std::vector<Triangle>& triangles, int maxLeafSize = 4)
// {
//     BVHNode* node = new BVHNode(triangles, maxLeafSize);

//     if (node->isLeaf) return node;

//     int axis = node->box.longest_axis();

//     std::vector<Triangle> sortedTriangles = triangles;
//     std::sort(sortedTriangles.begin(), sortedTriangles.end(), [axis](const Triangle& a, const Triangle& b) {
//         return a.centroid[axis] < b.centroid[axis];
//     });

//     size_t mid = sortedTriangles.size() / 2;
//     std::vector<Triangle> left(sortedTriangles.begin(), sortedTriangles.begin() + mid);
//     std::vector<Triangle> right(sortedTriangles.begin() + mid, sortedTriangles.end());

//     node->left  = buildBVH_CM(left, maxLeafSize);
//     node->right = buildBVH_CM(right, maxLeafSize);

//     node->triangles.clear();

//     return node;
// }

// BVHNode* buildBVH_SAH(const std::vector<Triangle>& triangles, int nBuckets=12)
// {   
//     BVHNode* node = new BVHNode(triangles,4);
    
//     if (node->isLeaf) return node;

//     int axis = node->box.longest_axis();
//     std::vector<Triangle> sortedTriangles = triangles;
//     std::sort(sortedTriangles.begin(), sortedTriangles.end(), [axis](const Triangle& a, const Triangle& b) {
//         return a.centroid[axis] < b.centroid[axis];
//     });

//     float cmin = sortedTriangles[0].centroid[axis];
//     float cmax = sortedTriangles.back().centroid[axis];
//     float extent = cmax - cmin;
//     if (extent < 1e-8f) extent = 1e-8f;
//     int best_split = 0;
//     if (extent==0){
//         best_split = sortedTriangles.size()/2;
//     }
//     else{
//         std::vector<Bucket> buckets(nBuckets);
//         for(Triangle t : sortedTriangles)
//         {
//             int b = nBuckets*(t.centroid[axis]-cmin) / extent;
//             if (b == nBuckets) b = nBuckets - 1;
//             buckets[b].count++;
//         }

//         int c = 0;
//         for(Bucket &b : buckets)
//         {   
//             b.box.update_box(sortedTriangles,c,c+b.count);
//             c +=b.count;
//         }

//         float minCost = INFINITY;
//         int count = 0;
//         for(int splits = 1;splits<buckets.size();++splits)
//         {
//             AABB left_box = buckets[0].box;
//             for (int i = 1;i<splits;++i)
//             {
//                 left_box = AABB(left_box,buckets[i].box);
//             }
//             AABB right_box = buckets[splits].box;
//             for (int i = splits+1;i<buckets.size();++i)
//             {
//                 right_box = AABB(right_box,buckets[i].box);
//             }
//             float cost = count*left_box.surface_area()+(sortedTriangles.size()-count)*right_box.surface_area();

//             count += buckets[splits-1].count;
//             if(cost < minCost){
//                 minCost = cost;
//                 best_split = count;
//             }
//         }
//     }

//     std::vector<Triangle> left(sortedTriangles.begin(), sortedTriangles.begin() + best_split);
//     std::vector<Triangle> right(sortedTriangles.begin() + best_split, sortedTriangles.end());

//     node->left  = buildBVH_SAH(left);
//     node->right = buildBVH_SAH(right);

//     return node;
// }

// bool intersectAABB(const Ray& ray, const AABB& box)
// {
//     float t_x0 = (box.l[0] - ray.point.x)/(ray.direction.x);
//     float t_x1 = (box.u[0] - ray.point.x)/(ray.direction.x);
//     if (t_x0 > t_x1) std::swap(t_x0, t_x1);

//     float tEnter = t_x0;
//     float tExit  = t_x1;

//     float t_y0 = (box.l[1] - ray.point.y)/(ray.direction.y);
//     float t_y1 = (box.u[1] - ray.point.y)/(ray.direction.y);
//     if (t_y0 > t_y1) std::swap(t_y0, t_y1);

//     tEnter = std::max(tEnter, t_y0);
//     tExit  = std::min(tExit,  t_y1);

//     if (tEnter > tExit) return false;

//     float t_z0 = (box.l[2] - ray.point.z)/(ray.direction.z);
//     float t_z1 = (box.u[2] - ray.point.z)/(ray.direction.z);
//     if (t_z0 > t_z1) std::swap(t_z0, t_z1);

//     tEnter = std::max(tEnter, t_z0);
//     tExit  = std::min(tExit,  t_z1);

//     if (tEnter > tExit) return false;

//     return tExit >= 0.0f;
// }

// float moeller_trumbore(const Ray& ray, const Triangle& tri, float& u, float& v)
// {
//     const float EPSILON = 1e-8f;
//     vec3 edge1 = tri.v[1] - tri.v[0];
//     vec3 edge2 = tri.v[2] - tri.v[0];

//     vec3 h = glm::cross(ray.direction, edge2);
//     float a = glm::dot(edge1, h);

//     if (a > -EPSILON && a < EPSILON)
//         return -1.0f; // Ray is parallel to the triangle

//     float f = 1.0f / a;
//     vec3 s = ray.point - tri.v[0];
//     u = f * glm::dot(s, h);

//     if (u < 0.0f || u > 1.0f)
//         return -1.0f;

//     vec3 q = glm::cross(s, edge1);
//     v = f * glm::dot(ray.direction, q);

//     if (v < 0.0f || u + v > 1.0f)
//         return -1.0f;

//     float t = f * glm::dot(edge2, q);

//     if (t > EPSILON) // ray intersection
//         return t;

//     return -1.0f;
// }

// std::vector<vec3> get_image_plane(float distance, Ray cam, uint16_t width,uint16_t height)
// {
//     vec3 center = cam.at(distance);
//     vec3 v1 = normalize(cross(cam.direction,{0.f,0.f,1.f}))*float(width)/100.0f;
//     vec3 v2 = normalize(cross(cam.direction,v1))*float(height)/100.0f;
//     center = center - 0.5f*v1 - 0.5f*v2;
//     return {center,v1,v2};
// }

// HitRecord traverseBVH(const BVHNode* root, const Ray& ray, bool shadow_ray)
// {
//     HitRecord best_hit(INFINITY, nullptr);
//     float u,v = -1.f;
//     std::stack<const BVHNode*> stack;
//     stack.push(root);
//     LOG(heatMap.intersects = 0;);

//     while (!stack.empty()) {
//         const BVHNode* node = stack.top();
//         stack.pop();
//         LOG(heatMap.intersects += 1;);

//         if (!intersectAABB(ray, node->box))
//             continue;

//         if (node->isLeaf) {
//             for (const Triangle& t : node->triangles){
//                 float t_hit = moeller_trumbore(ray, t,u,v);
//                 if (t_hit < best_hit.t && t_hit >= 0.f) 
//                 {
//                     best_hit.triangle = &t;
//                     best_hit.t = t_hit;
//                     best_hit.u = u;
//                     best_hit.v = v;
//                     if (shadow_ray) return best_hit;
//                 }
//             }
//         } else {
//             stack.push(node->left);
//             stack.push(node->right);
//         }
//     }
//     return best_hit;
// }

inline float reinhard(float x)
{
    return x / (1.0f + x);
}

inline float aces(float v)
{
    float a = v * (v + 0.0245786) - 0.000090537;
    float b = v * (0.983729 * v + 0.4329510) + 0.238081;
    return a/b;
}

std::vector<uint8> tone_mapping(const std::vector<float>& float_pixels, float yamma)
{
    std::vector<uint8> pixels((float_pixels.size()*4)/3);
    for(size_t idx=0;idx<pixels.size()/4;idx++)
    {
        for(int c=0;c<3;++c){
            uint8 v = uint8(std::pow(reinhard(float_pixels[idx*3 + c]),yamma)*255); 
            pixels[idx*4 + c] = v;
        }
        pixels[idx*4 + 3] = 255;
    }
    return pixels;
}

int main()
{ 
    // Ray cam({30.f,40.f,10.f},{0.f,-15.f,5.f});
    // float distance = 5.f;
    // uint16_t width = 1080;
    // uint16_t height = 720;
    // std::vector<vec3> image_plane = get_image_plane(distance,cam,width,height);
    // vec3 img_point = image_plane[0];
    // vec3 veci = image_plane[1];
    // vec3 vecj = image_plane[2];
    // std::vector<float> float_pixels(width * height * 3, 0);
    // heatMap.heatpixels.reserve(float_pixels.size()/3);

    // LightSource point_light({2000.f,1000.f,1000.f},{1.f,1.0f,1.0f});
    // vec3 ambient_light = {1.0,1.0,1.0};
    // float alpha = 4;
    // float beta = 1;
    // float gamma = 0.001;
    // float m = 10;

    //NEW MAIN


    std::vector<Model> models;
    std::vector<std::string> models_to_load = {"plane", "trees9"};
    for(std::string m : models_to_load) models.emplace_back(m,"/home/lukas/simple-raytracer/scene/");
    std::cout << "All Models loaded! Scene is getting set up..." << std::endl;

    Scene scene(models,"/home/lukas/simple-raytracer/scene/setup.yaml");
    scene.update_data();
    std::cout << "Scene setup successful! (contains " << scene.triangles.size() << " triangles) Starting BVH Building..." << std::endl;

    BVH bvh(scene.triangles);
    std::cout << "Building is finished! Starting Rendering..." << std::endl;

    //END
    
    auto start = std::chrono::high_resolution_clock::now();


    int height = scene.camera.resolution.x;
    int width = scene.camera.resolution.y;
    std::vector<float> float_pixels(width * height * 3, 0);

    for(int j=0;j<height;++j)
    {
        for(int i=0;i<width;++i)
        {
            vec3 pos = scene.camera.image_plane.origin + scene.camera.image_plane.w*float(i)/float(width) + scene.camera.image_plane.h*float(j)/float(height);
            Ray ray(scene.camera.position, pos);
            unsigned idx = (j * width + i) * 3;
            
            HitRecord hit = bvh.traverse_BVH(ray, false);
            LOG(heatMap.heatpixels.push_back(heatMap.intersects));
            if(hit.t!=INFINITY)
            {
                Ray shadow_ray(ray.at(hit.t)+hit.triangle->normal*0.0005f,scene.lights.front().position);
                HitRecord shadow_hit = bvh.traverse_BVH(shadow_ray, true);

                float intensity = 0.f;
                vec3 v(0.f);
                vec3 color_diffuse(0.f);
                if(shadow_hit.t==INFINITY) {
                    intensity = std::max(0.0f,dot(shadow_ray.direction, hit.triangle->normal));
                    v = shadow_ray.direction - 2*dot(hit.triangle->normal,-shadow_ray.direction)*hit.triangle->normal;
                }

                int tex_id = scene.materials[hit.triangle->mat_id].diffuseTex;
                if(tex_id >=0)
                {
                    vec2 uv_values = hit.u*hit.triangle->uv[0] + hit.v*hit.triangle->uv[1] + (1-hit.u-hit.v) * hit.triangle->uv[2];
                    Texture tex = scene.textures[tex_id];

                    int x = clamp(int(uv_values.x * (tex.width  - 1)), 0, tex.width  - 1);
                    int y = clamp(int(uv_values.y * (tex.height - 1)), 0, tex.height - 1);
                    int tex_idx = (y * tex.width + x) * tex.channels;

                    color_diffuse = {tex.data[tex_idx + 0],tex.data[tex_idx + 1],tex.data[tex_idx + 2]};
                }
                else
                {
                    color_diffuse = scene.materials[hit.triangle->mat_id].diffuseColor;
                }
                
                for(int c=0;c<3;++c){
                    float_pixels[idx + c] = color_diffuse[c]*(alpha*point_light.color[c]*intensity+beta*ambient_light[c])
                    +(gamma*point_light.color[c]*intensity*std::pow(dot(-ray.direction,v),m));                   
                }
            }
        }
    }

    auto pixels = tone_mapping(float_pixels,1.6);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Total render time: " << elapsed.count() << " seconds\n";

    auto window = sf::RenderWindow(sf::VideoMode({width, height}), "RaytracerPOG");
    window.setFramerateLimit(60);
    

    sf::Texture texture(sf::Vector2u(width, height));
    LOG(texture.update(generateHeatmap(heatMap.heatpixels).data()););
    if (HEATMAP==false) texture.update(pixels.data());

    sf::Sprite sprite(texture);


    while (window.isOpen())
    {
        while (const std::optional event = window.pollEvent())
        {
            if (event->is<sf::Event::Closed>())
            {
                window.close();
            }
        }
        
        window.clear();
        window.draw(sprite);
        window.display();
    }
}
