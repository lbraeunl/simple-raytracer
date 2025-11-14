#pragma once
#include <glm/glm.hpp>
#include <iostream>

struct Ray
{
    glm::vec3 point;
    glm::vec3 direction;

    Ray();
    Ray(const glm::vec3& p1, const glm::vec3& p2);

    glm::vec3 at(float t) const;
};

struct Triangle
{
    glm::vec3 a, b, c;
    glm::vec3 color;
    glm::vec3 centroid;
    glm::vec3 normal;

    Triangle();
    Triangle(const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& v3, const glm::vec3& color);
};

struct AABB
{
    glm::vec3 l;
    glm::vec3 u;
    AABB();
    AABB(const std::vector<Triangle>& triangles);

    uint8_t longest_axis() const;
    float surface_area() const;
};

struct BVHNode {
    AABB box;
    BVHNode* left;
    BVHNode* right;
    bool isLeaf;
    std::vector<Triangle> triangles;

    BVHNode(const std::vector<Triangle>& triangles,int max_leaf_size);
};

struct HitRecord {
    float t;
    const Triangle* triangle;

    HitRecord(const float& t, Triangle* triangle);
};

struct LightSource {
    glm::vec3 position;
    glm::vec3 color;

    LightSource(const glm::vec3& position, const glm::vec3& color);
};