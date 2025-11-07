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

    Triangle();
    Triangle(const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& v3, const glm::vec3& color);

    glm::vec3 normal() const;
};

struct AABB
{
    float x_l, x_u, y_l, y_u, z_l, z_u;
    //glm::vec3 centroid;
    AABB();
    AABB(const std::vector<Triangle>& triangles);

    uint8_t longest_axis() const;
};

struct BVHNode {
    AABB box;
    BVHNode* left;
    BVHNode* right;
    bool isLeaf;
    std::vector<Triangle> triangles;

    BVHNode(const std::vector<Triangle>& triangles,int max_leaf_size);
};
