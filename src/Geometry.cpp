#define GLM_ENABLE_EXPERIMENTAL
#include "Geometry.hpp"
#include <glm/gtx/norm.hpp>
#include <glm/ext/scalar_common.hpp>
#include "stb_image.h"
#include "logger.hpp"

//Miscellaneous

uint8_t argmax(glm::vec3 v)
{
    if (v.x>= v.y && v.x >= v.z) return (uint8_t)0;
    if (v.y >=v.z)return (uint8_t)1;
    return (uint8_t)2;
}

// Ray

Ray::Ray() : point(0.0f), direction(0.0f, 1.0f, 0.0f) {}

Ray::Ray(const glm::vec3& p1, const glm::vec3& p2)
    : point(p1), direction(glm::normalize(p2 - p1)) {}

glm::vec3 Ray::at(float t) const {
    return point + t * direction;
}


// Triangle

Triangle::Triangle()
    : centroid(0.f), normal(0.f), mat_id(-1) {}

Triangle::Triangle(glm::vec3 vertices[3], glm::vec2 uvs[3], glm::vec3 normals[3], int mat_id)
{
    for(int i=0;i<3;++i)
    {
        v[i] = vertices[i];
        uv[i] = uvs[i];
        n[i] = normals[i];
    }
    update();
    mat_id = mat_id;
}

Triangle::Triangle(glm::vec3 a,glm::vec3 b,glm::vec3 c,int mat_id)
{
    v[0]=a;
    v[1]=b;
    v[2]=c;
    update();
    mat_id = mat_id;
}

void Triangle::update()
{
    centroid = (v[0] + v[1] + v[2]) / 3.0f;
    normal = glm::normalize(glm::cross(v[1] - v[0], v[2] - v[0]));
    glm::vec3 e1 = v[1] - v[0];
    glm::vec3 e2 = v[2] - v[0];
    area = 0.5f * glm::length(glm::cross(e1, e2));
}

void Triangle::print() const
{
    std::cout << "vertex 1:" << v[0].x << "," << v[0].y << "," << v[0].z << std::endl;
    std::cout << "vertex 2:" << v[1].x << "," << v[1].y << "," << v[1].z << std::endl;
    std::cout << "vertex 3:" << v[2].x << "," << v[2].y << "," << v[2].z << std::endl;
    std::cout << "normal:" << normal.x << "," << normal.y << "," << normal.z << std::endl;
}

glm::vec3 Triangle::interpolate_normal(float u, float v) const
    {
        if (hasVertexNormals) 
        {
            return glm::normalize(u * n[1] + v * n[2] + (1.0f - u - v) * n[0]);
        }
        return normal;
    }

bool Triangle::triangle_intersect(const Ray& ray, HitRecord& hit) const
{
    LOG(heat_map.intersects+=1;);
    const float EPSILON = 1e-8f;
    glm::vec3 edge1 = v[1] - v[0];
    glm::vec3 edge2 = v[2] - v[0];

    glm::vec3 h = glm::cross(ray.direction, edge2);
    float a = glm::dot(edge1, h);

    if (/*a > -EPSILON && */ a < EPSILON)
        return false;

    float f = 1.0f / a;
    glm::vec3 s = ray.point - v[0];
    hit.u = f * glm::dot(s, h);

    if (hit.u < 0.0f || hit.u > 1.0f)
        return false;

    glm::vec3 q = glm::cross(s, edge1);
    hit.v = f * glm::dot(ray.direction, q);

    if (hit.v < 0.0f || hit.u + hit.v > 1.0f)
        return false;

    hit.t = f * glm::dot(edge2, q);

    if (hit.t <= EPSILON)
        return false;

    hit.triangle = this;
    hit.position = ray.at(hit.t);
    hit.normal = interpolate_normal(hit.u, hit.v);
    return -1.0f;
}


// AABB

AABB::AABB()
    : l(INFINITY,INFINITY,INFINITY), u(-INFINITY,-INFINITY,-INFINITY) {}

// AABB::AABB(const std::vector<Triangle>& triangles)
//     {
//         update_box(triangles,0,triangles.size());
//     };

AABB::AABB(const Triangle& t)
{
    l = glm::min(t.v[0],t.v[1],t.v[2]);
    u = glm::max(t.v[0],t.v[1],t.v[2]);
}

AABB::AABB(const AABB b1, const AABB b2)
{
    l = glm::min(b1.l,b2.l);
    u = glm::max(b1.u,b2.u);
}

float AABB::surface_area() const
{
    return (u[0]-l[0])*(u[1]-l[1])+(u[0]-l[0])*(u[2]-l[2])+(u[1]-l[1])*(u[2]-l[2]);
}

// void AABB::update_box(const std::vector<Triangle>& sorted_triangles,int start,int end)
// {
//     l = {INFINITY,INFINITY,INFINITY};
//     u = {-INFINITY,-INFINITY,-INFINITY};
//     for (int i = start;i<end;++i)
//     {         
//         expand(sorted_triangles[i]);
//     }  
// }

void AABB::expand(AABB b)
{
    for(int i=0;i<3;++i)
    {
        if(b.l[i]<l[i]) l[i] = b.l[i];
        if(b.u[i]>u[i]) u[i] = b.u[i];
    }
}

uint8_t AABB::longest_axis() const
{
    return argmax(glm::vec3(u-l));
}

bool AABB::box_intersect(const Ray& ray, float t_max) const
{
    LOG(heat_map.intersects+=1;);
    float t_x0 = (l[0] - ray.point.x)/(ray.direction.x);
    float t_x1 = (u[0] - ray.point.x)/(ray.direction.x);
    if (t_x0 > t_x1) std::swap(t_x0, t_x1);

    float tEnter = t_x0;
    float tExit  = t_x1;

    float t_y0 = (l[1] - ray.point.y)/(ray.direction.y);
    float t_y1 = (u[1] - ray.point.y)/(ray.direction.y);
    if (t_y0 > t_y1) std::swap(t_y0, t_y1);

    tEnter = std::max(tEnter, t_y0);
    tExit  = std::min(tExit,  t_y1);

    if (tEnter > tExit) return false;

    float t_z0 = (l[2] - ray.point.z)/(ray.direction.z);
    float t_z1 = (u[2] - ray.point.z)/(ray.direction.z);
    if (t_z0 > t_z1) std::swap(t_z0, t_z1);

    tEnter = std::max(tEnter, t_z0);
    tExit  = std::min(tExit,  t_z1);

    if (tEnter > tExit) return false;

    if (tEnter > t_max) return false;

    return tExit >= 0.0f;
}


// HitRecord

HitRecord::HitRecord() {}

HitRecord::HitRecord(const float& t, Triangle* triangle):t(t),triangle(triangle), u(-1.f), v(-1.f) {}


//Rectangle

Rectangle::Rectangle(const glm::vec3 origin, const glm::vec3 w, const glm::vec3 h): origin(origin), w(w), h(h) {}

