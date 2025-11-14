#include <SFML/Graphics.hpp>
#include "tiny_obj_loader.h"
#include <iostream>
#include <glm/vec3.hpp>
#include "Geometry.hpp"
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/norm.hpp>
#include <algorithm>
#include <stack>
#include <chrono>
#include <tinycolormap.hpp>

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

BVHNode* buildBVH_CM(const std::vector<Triangle>& triangles, int maxLeafSize = 4)
{
    BVHNode* node = new BVHNode(triangles, maxLeafSize);

    if (node->isLeaf) return node;

    int axis = node->box.longest_axis();

    std::vector<Triangle> sortedTriangles = triangles;
    std::sort(sortedTriangles.begin(), sortedTriangles.end(), [axis](const Triangle& a, const Triangle& b) {
        return a.centroid[axis] < b.centroid[axis];
    });

    size_t mid = sortedTriangles.size() / 2;
    std::vector<Triangle> left(sortedTriangles.begin(), sortedTriangles.begin() + mid);
    std::vector<Triangle> right(sortedTriangles.begin() + mid, sortedTriangles.end());

    node->left  = buildBVH_CM(left, maxLeafSize);
    node->right = buildBVH_CM(right, maxLeafSize);

    node->triangles.clear();
    node->isLeaf = false;

    return node;
}

// BVHNode* buildBVH_SAH(const std::vector<Triangle>& triangles, int nBuckets=12)
// {
//     BVHNode* node = new BVHNode(triangles, maxLeafSize);

//     if (node->isLeaf) return node;

//     int axis = node->box.longest_axis();

//     std::vector<Triangle> sortedTriangles = triangles;
//     std::sort(sortedTriangles.begin(), sortedTriangles.end(), [axis](const Triangle& a, const Triangle& b) {
//         return a.centroid[axis] < b.centroid[axis];
//     });

//     std::vector<BVHNode> buckets(nBuckets);

//     for(Triangle t : node->triangles)
//     {
//         node->box
//            t.centroid[axis]
//     }

    



// }


bool intersectAABB(const Ray& ray, const AABB& box)
{
    
    
    
    
    float t_x0 = (box.l[0] - ray.point.x)/(ray.direction.x);
    float t_x1 = (box.u[0] - ray.point.x)/(ray.direction.x);
    if (t_x0 > t_x1) std::swap(t_x0, t_x1);

    float tEnter = t_x0;
    float tExit  = t_x1;

    float t_y0 = (box.l[1] - ray.point.y)/(ray.direction.y);
    float t_y1 = (box.u[1] - ray.point.y)/(ray.direction.y);
    if (t_y0 > t_y1) std::swap(t_y0, t_y1);

    tEnter = std::max(tEnter, t_y0);
    tExit  = std::min(tExit,  t_y1);

    if (tEnter > tExit) return false;

    float t_z0 = (box.l[2] - ray.point.z)/(ray.direction.z);
    float t_z1 = (box.u[2] - ray.point.z)/(ray.direction.z);
    if (t_z0 > t_z1) std::swap(t_z0, t_z1);

    tEnter = std::max(tEnter, t_z0);
    tExit  = std::min(tExit,  t_z1);

    if (tEnter > tExit) return false;

    return tExit >= 0.0f;
}

float moeller_trumbore(const Ray& ray, const Triangle& tri)
{
    const float EPSILON = 1e-8f;
    vec3 edge1 = tri.b - tri.a;
    vec3 edge2 = tri.c - tri.a;

    vec3 h = glm::cross(ray.direction, edge2);
    float a = glm::dot(edge1, h);

    if (a > -EPSILON && a < EPSILON)
        return -1.0f; // Ray is parallel to the triangle

    float f = 1.0f / a;
    vec3 s = ray.point - tri.a;
    float u = f * glm::dot(s, h);

    if (u < 0.0f || u > 1.0f)
        return -1.0f;

    vec3 q = glm::cross(s, edge1);
    float v = f * glm::dot(ray.direction, q);

    if (v < 0.0f || u + v > 1.0f)
        return -1.0f;

    float t = f * glm::dot(edge2, q);

    if (t > EPSILON) // ray intersection
        return t;

    return -1.0f;
}

std::vector<Triangle> load_object(std::string filename, std::string directory="/home/lukas/simple-raytracer")
{
    
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;

	std::string warn;
	std::string err;
    
    tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filename.c_str(),directory.c_str());

    if (!warn.empty()) {
		std::cout << "WARN: " << warn << std::endl;
	}
    
	if (!err.empty()) {
		std::cerr << err << std::endl;
		return {};
	}

    std::vector<Triangle> triangles;
    size_t total_faces = 0;

    for (const auto& shape : shapes) {
        total_faces += shape.mesh.num_face_vertices.size();
    }
    triangles.reserve(total_faces);

    for (const auto& shape : shapes) 
    {
        for (size_t f = 0; f < shape.mesh.num_face_vertices.size(); f++) 
        {
            if (shape.mesh.num_face_vertices[f] != 3) {
                std::cout << "Warning:" << filename << " contains a non-triangular shape that is skipped." << std::endl;
                continue;
            }

            tinyobj::index_t idx0 = shape.mesh.indices[3*f+0];
            tinyobj::index_t idx1 = shape.mesh.indices[3*f+1];
            tinyobj::index_t idx2 = shape.mesh.indices[3*f+2];
            int mat_id = shape.mesh.material_ids[f];


            vec3 v0(attrib.vertices[3*idx0.vertex_index+0], attrib.vertices[3*idx0.vertex_index+1], attrib.vertices[3*idx0.vertex_index+2]);
            vec3 v1(attrib.vertices[3*idx1.vertex_index+0], attrib.vertices[3*idx1.vertex_index+1], attrib.vertices[3*idx1.vertex_index+2]);
            vec3 v2(attrib.vertices[3*idx2.vertex_index+0], attrib.vertices[3*idx2.vertex_index+1], attrib.vertices[3*idx2.vertex_index+2]);


            vec3 color(
                float(materials[mat_id].diffuse[0]),
                float(materials[mat_id].diffuse[1]),
                float(materials[mat_id].diffuse[2])
                );

            triangles.emplace_back(v0, v1, v2, color);
        }
    }

    return triangles;
}

std::vector<vec3> get_image_plane(float distance, Ray cam, uint16_t width,uint16_t height)
{
    vec3 center = cam.at(distance);
    vec3 v1 = normalize(cross(cam.direction,{0.f,0.f,1.f}))*float(width)/100.0f;
    vec3 v2 = normalize(cross(cam.direction,v1))*float(height)/100.0f;
    center = center - 0.5f*v1 - 0.5f*v2;
    return {center,v1,v2};
}

void addFloor(std::vector<Triangle>& triangles, uint16 size=400)
{
    triangles.push_back(Triangle({size*0.5f,size*(-0.5f),0.f},{size*0.5f, size*0.5f,0.f},{size*-0.5f, size*-0.5f,0.f},{0.4f,0.4f,0.4f}));
    triangles.push_back(Triangle({size*0.5f, size*0.5f,0.f},{size*(-0.5f),size*0.5f,0.f},{size*-0.5f, size*-0.5f,0.f},{0.4f,0.4f,0.4f}));
}

HitRecord traverseBVH(const BVHNode* root, const Ray& ray, bool shadow_ray)
{
    HitRecord best_hit(INFINITY, nullptr);
    std::stack<const BVHNode*> stack;
    stack.push(root);
    LOG(heatMap.intersects = 0;);

    while (!stack.empty()) {
        const BVHNode* node = stack.top();
        stack.pop();
        LOG(heatMap.intersects += 1;);

        if (!intersectAABB(ray, node->box))
            continue;

        if (node->isLeaf) {
            for (const Triangle& t : node->triangles){
                float t_hit = moeller_trumbore(ray, t);
                if (t_hit < best_hit.t && t_hit >= 0.f) 
                {
                    best_hit.triangle = &t;
                    best_hit.t = t_hit;
                    if (shadow_ray) return best_hit;
                }      
            }
        } else {
            stack.push(node->left);
            stack.push(node->right);
        }
    }
    return best_hit;
}

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
    Ray cam({100.f,40.f,40.f},{0.f,0.f,10.f});
    float distance = 5.f;
    uint16_t width = 1080;
    uint16_t height = 720;
    std::vector<vec3> image_plane = get_image_plane(distance,cam,width,height);
    vec3 img_point = image_plane[0];
    vec3 veci = image_plane[1];
    vec3 vecj = image_plane[2];
    std::vector<float> float_pixels(width * height * 3, 0);
    heatMap.heatpixels.reserve(float_pixels.size()/3);

    LightSource point_light({-2000.f,-1000.f,3000.f},{1.f,1.0f,1.0f});
    vec3 ambient_light = {1.0,1.0,1.0};
    float alpha = 4;
    float beta = 0.5;
    float gamma = 0.1;
    float m = 10;

    std::vector<Triangle> triangles = load_object("/home/lukas/simple-raytracer/tinker.obj");
    addFloor(triangles);

    BVHNode* node = buildBVH_CM(triangles);
    std::cout << "Building is finished! Starting Rendering..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    for(int j=0;j<height;++j)
    {
        for(int i=0;i<width;++i) 
        {
            vec3 pos = img_point+(veci*float(i))/float(width)+(vecj*float(j))/float(height);
            Ray ray(cam.point, pos);
            unsigned idx = (j * width + i) * 3;
            
            HitRecord hit = traverseBVH(node, ray, false);
            LOG(heatMap.heatpixels.push_back(heatMap.intersects));
            if(hit.t!=INFINITY)
            {
                Ray shadow_ray(ray.at(hit.t)+hit.triangle->normal*0.0005f,point_light.position);
                HitRecord shadow_hit = traverseBVH(node, shadow_ray, true);

                float intensity = 0.f;
                vec3 v(0.f);
                if(shadow_hit.t==INFINITY) {
                    intensity = std::max(0.0f,dot(shadow_ray.direction, hit.triangle->normal));
                    v = shadow_ray.direction - 2*dot(hit.triangle->normal,-shadow_ray.direction)*hit.triangle->normal;
                }

                for(int c=0;c<3;++c){
                    float_pixels[idx + c] = hit.triangle->color[c]*(alpha*point_light.color[c]*intensity+beta*ambient_light[c])
                    +(gamma*point_light.color[c]*intensity*std::pow(dot(-ray.direction,v),m));                   
                }
            }
        }
    }

    auto pixels = tone_mapping(float_pixels,1.7);

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
