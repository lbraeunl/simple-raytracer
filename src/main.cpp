#include <SFML/Graphics.hpp>
#include "tiny_obj_loader.h"
#include <iostream>
#include <glm/vec3.hpp>
#include "Geometry.hpp"
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/norm.hpp>
#include <stack>

using namespace glm;

BVHNode* buildBVH(const std::vector<Triangle>& triangles, int maxLeafSize = 4)
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

    node->left  = buildBVH(left, maxLeafSize);
    node->right = buildBVH(right, maxLeafSize);

    node->triangles.clear();
    node->isLeaf = false;

    return node;
}

bool intersectAABB(const Ray& ray, const AABB& box)
{
    float t_x0 = (box.x_l - ray.point.x)/(ray.direction.x);
    float t_x1 = (box.x_u - ray.point.x)/(ray.direction.x);
    if (t_x0 > t_x1) std::swap(t_x0, t_x1);

    float tEnter = t_x0;
    float tExit  = t_x1;

    float t_y0 = (box.y_l - ray.point.y)/(ray.direction.y);
    float t_y1 = (box.y_u - ray.point.y)/(ray.direction.y);
    if (t_y0 > t_y1) std::swap(t_y0, t_y1);

    tEnter = std::max(tEnter, t_y0);
    tExit  = std::min(tExit,  t_y1);

    if (tEnter > tExit) return false;

    float t_z0 = (box.z_l - ray.point.z)/(ray.direction.z);
    float t_z1 = (box.z_u - ray.point.z)/(ray.direction.z);
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
                uint8_t(materials[mat_id].diffuse[0]*255),
                uint8_t(materials[mat_id].diffuse[1]*255),
                uint8_t(materials[mat_id].diffuse[2]*255)
                );

            triangles.emplace_back(v0, v1, v2, color);
        }
    }

    return triangles;
}

std::vector<vec3> get_image_plane(float distance, Ray cam, float size=10.f)
{
    vec3 center = cam.at(distance);
    vec3 v1 = normalize(cross(cam.direction,{0.f,0.f,1.f}))*size;
    vec3 v2 = normalize(cross(cam.direction,v1))*size;
    center = center - 0.5f*v1 - 0.5f*v2;
    return {center,v1,v2};
}

HitRecord traverseBVH(const BVHNode* root, const Ray& ray)
{
    HitRecord best_hit(INFINITY, vec3(0.f));
    std::stack<const BVHNode*> stack;
    stack.push(root);

    while (!stack.empty()) {
        const BVHNode* node = stack.top();
        stack.pop();

        if (!intersectAABB(ray, node->box))
            continue;

        if (node->isLeaf) {
            for (const Triangle& t : node->triangles){
                float t_hit = moeller_trumbore(ray, t);
                if (t_hit < best_hit.t && t_hit >= 0.f) 
                {
                    best_hit.color = t.color;
                    best_hit.t = t_hit;
                }      
            }
        } else {
            stack.push(node->left);
            stack.push(node->right);
        }
    }
    return best_hit;
}

int main()
{ 
    Ray cam({-100.f,-40.f,40.f},{0.f,0.f,10.f});
    float distance = 5.f;
    
    std::vector<vec3> image_plane = get_image_plane(distance,cam);
    vec3 img_point = image_plane[0];
    vec3 veci = image_plane[1];
    vec3 vecj = image_plane[2];

    uint16_t width = 1080;
    uint16_t height = 1080;
    std::vector<uint8_t> pixels(width * height * 4, 0);

    std::vector<Triangle> triangles = load_object("/home/lukas/simple-raytracer/tinker.obj");

    BVHNode* node = buildBVH(triangles);
    std::cout << "Building is finished! Starting Rendering..." << std::endl;

    for(int i=0;i<width;++i)
    {
        std::cout <<"\r"<< (float(i)/width)*100 << "%";
        for(int j=0;j<height;++j) 
        {
            vec3 pos = img_point+(veci*float(i))/float(width)+(vecj*float(j))/float(height);
            Ray ray(cam.point, pos);
            
            HitRecord hit = traverseBVH(node, ray);

            unsigned idx = (j * width + i) * 4;
            pixels[idx + 0] = hit.color[0];
            pixels[idx + 1] = hit.color[1];
            pixels[idx + 2] = hit.color[2];
            pixels[idx + 3] = 255;
        }
    }

    
    auto window = sf::RenderWindow(sf::VideoMode({width, height}), "RaytracerPOG");
    window.setFramerateLimit(60);
    

    sf::Texture texture(sf::Vector2u(width, height));
    texture.update(pixels.data());

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
