#define GLM_ENABLE_EXPERIMENTAL
#include "Geometry.hpp"
#include <glm/gtx/norm.hpp>
#include <glm/ext/scalar_common.hpp>

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
    : a(0.0f), b(1.0f, 0.0f, 0.0f), c(0.0f, 1.0f, 0.0f), color(0,0,0), centroid(0,0,0), normal(0,0,0) {}

Triangle::Triangle(const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& v3, const glm::vec3& color)
    : a(v1), b(v2), c(v3), color(color) {
        centroid = (a + b + c) / 3.0f;
        normal = glm::normalize(glm::cross(b - a, c - a));
    }


// AABB

AABB::AABB()
    : l(0.0f,0.f,0.f), u(0.f,0.f,0.f) {}

float AABB::surface_area() const
{
    return (u[0]-l[0])*(u[1]-l[1])+(u[0]-l[0])*(u[2]-l[2])+(u[1]-l[1])*(u[2]-l[2]);
}

AABB::AABB(const std::vector<Triangle>& triangles)
    {
        l = {INFINITY,INFINITY,INFINITY};
        u = {-INFINITY,-INFINITY,-INFINITY};
        for (const Triangle& t :triangles)
        {
            glm::vec3 min = glm::min(t.a,t.b,t.c);
            glm::vec3 max = glm::max(t.a,t.b,t.c);
            for(int i=0;i<3;++i) 
            {
                if(min[i]<l[i]) l[i] = min[i];
                if(max[i]>u[i]) u[i] = max[i];
            }
        }
    };

uint8_t AABB::longest_axis() const
{
    return argmax(glm::vec3(u-l));
}


// BVHNode

BVHNode::BVHNode(const std::vector<Triangle>& triangles, int max_leaf_size):triangles(triangles), box(AABB(triangles)) {
    isLeaf = triangles.size()<= max_leaf_size;
}


// HitRecord

HitRecord::HitRecord(const float& t, Triangle* triangle):t(t),triangle(triangle){}


//LightSource

LightSource::LightSource(const glm::vec3& position, const glm::vec3& color): position(position), color(color){}