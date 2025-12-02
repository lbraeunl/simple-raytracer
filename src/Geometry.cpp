#define GLM_ENABLE_EXPERIMENTAL
#include "Geometry.hpp"
#include <glm/gtx/norm.hpp>
#include <glm/ext/scalar_common.hpp>
#include "stb_image.h"

//Miscellaneous

uint8_t argmax(glm::vec3 v)
{
    if (v.x>= v.y && v.x >= v.z) return (uint8_t)0;
    if (v.y >=v.z)return (uint8_t)1;
    return (uint8_t)2;
}

// Line

Ray::Ray() : point(0.0f), direction(0.0f, 1.0f, 0.0f) {}

Ray::Ray(const glm::vec3& p1, const glm::vec3& p2)
    : point(p1), direction(glm::normalize(p2 - p1)) {}

glm::vec3 Ray::at(float t) const {
    return point + t * direction;
}


// Triangle

Triangle::Triangle()
    : centroid(0.f), normal(0.f), mat_id(-1) {}

Triangle::Triangle(glm::vec3 vertices[3], glm::vec2 uvs[3], int mat_id)
{
    for(int i=0;i<3;++i)
    {
        v[i] = vertices[i];
        uv[i] = uvs[i];
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
}

void Triangle::yz_swap()
{   
    for(int i=0;i<3;++i)
    {
        std::swap(v[i][1],v[i][2]);
    }
}

void Triangle::print() const
{
    std::cout << "vertice1:" << v[0].x << "," << v[0].y << "," << v[0].z << std::endl;
    std::cout << "vertice2:" << v[1].x << "," << v[1].y << "," << v[1].z << std::endl;
    std::cout << "vertice3:" << v[2].x << "," << v[2].y << "," << v[2].z << std::endl;
    std::cout << "normal:" << normal.x << "," << normal.y << "," << normal.z << std::endl;
}


// AABB

AABB::AABB()
    : l(0.0f,0.f,0.f), u(0.f,0.f,0.f) {}

AABB::AABB(const std::vector<Triangle>& triangles)
    {
        update_box(triangles,0,triangles.size());
    };

AABB::AABB(const AABB b1, const AABB b2)
{
    l = glm::min(b1.l,b2.l);
    u = glm::max(b1.u,b2.u);
}

float AABB::surface_area() const
{
    return (u[0]-l[0])*(u[1]-l[1])+(u[0]-l[0])*(u[2]-l[2])+(u[1]-l[1])*(u[2]-l[2]);
}

void AABB::update_box(const std::vector<Triangle>& sorted_triangles,int start,int end)
{
    l = {INFINITY,INFINITY,INFINITY};
    u = {-INFINITY,-INFINITY,-INFINITY};
    for (int i = start;i<end;++i)
    {   
        const Triangle& t = sorted_triangles[i];
        glm::vec3 min = glm::min(t.v[0],t.v[1],t.v[2]);
        glm::vec3 max = glm::max(t.v[0],t.v[1],t.v[2]);
        for(int i=0;i<3;++i) 
        {
            if(min[i]<l[i]) l[i] = min[i];
            if(max[i]>u[i]) u[i] = max[i];
        }
    }  
}

uint8_t AABB::longest_axis() const
{
    return argmax(glm::vec3(u-l));
}


// BVHNode

BVHNode::BVHNode(const std::vector<Triangle>& triangles, int max_leaf_size):triangles(triangles), box(AABB(triangles)) {
    isLeaf = triangles.size()<= max_leaf_size;
}

BVHNode::BVHNode(const std::vector<Triangle>& triangles,AABB box):box(box),triangles(triangles) {}

BVHNode::BVHNode(): box(AABB()) {
}


// HitRecord

HitRecord::HitRecord(const float& t, Triangle* triangle):t(t),triangle(triangle), u(-1.f), v(-1.f) {}


// Image Plane
