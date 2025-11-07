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
    : a(0.0f), b(1.0f, 0.0f, 0.0f), c(0.0f, 1.0f, 0.0f), color(0,0,0), centroid(0,0,0) {}

Triangle::Triangle(const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& v3, const glm::vec3& color)
    : a(v1), b(v2), c(v3), color(color) {
        centroid = (a + b + c) / 3.0f;
    }

glm::vec3 Triangle::normal() const {
    return glm::normalize(glm::cross(b - a, c - a));
}


// AABB

AABB::AABB()
    : x_l(0.0f), x_u(0.0f), y_l(0.0f), y_u(0.0f), z_l(0.0f),z_u(0.0f) {}

AABB::AABB(const std::vector<Triangle>& triangles)
    {
        x_l = INFINITY;
        y_l = INFINITY;
        z_l = INFINITY;
        x_u = -INFINITY;
        y_u = -INFINITY;
        z_u = -INFINITY;
        for (const Triangle& t :triangles)
        {
            glm::vec3 min = glm::min(t.a,t.b,t.c);
            if (min.x<x_l) x_l = min.x;
            if (min.y<y_l) y_l = min.y;
            if (min.z<z_l) z_l = min.z;
            glm::vec3 max = glm::max(t.a,t.b,t.c);
            if (max.x>x_u) x_u = max.x;
            if (max.y>y_u) y_u = max.y;
            if (max.z>z_u) z_u = max.z;
        }
    };

uint8_t AABB::longest_axis() const
{
    return argmax(glm::vec3(x_u-x_l,y_u-y_l,z_u-z_l));
}


// BVHNode

BVHNode::BVHNode(const std::vector<Triangle>& triangles, int max_leaf_size):triangles(triangles), box(AABB(triangles)) {
    isLeaf = triangles.size()<= max_leaf_size;
}


// HitRecord

HitRecord::HitRecord(const float& t, glm::vec3 color):t(t),color(color){}