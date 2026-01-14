#pragma once
#include <glm/glm.hpp>
#include <iostream>

struct Triangle;

struct HitRecord {
    float t;
    float u,v;
    const Triangle* triangle;  //refactor to only contain information and not whole triangle

    HitRecord();
    HitRecord(const float& t, Triangle* triangle);
};

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
    glm::vec3 v[3];
    glm::vec2 uv[3];
    glm::vec3 n[3];
    glm::vec3 centroid;
    glm::vec3 normal;
    int mat_id;

    Triangle();
    Triangle(glm::vec3 a,glm::vec3 b,glm::vec3 c,int mat_id);
    Triangle(glm::vec3 vertices[3], glm::vec2 uvs[3], glm::vec3 normals[3], int mat_id);

    void update();
    void print() const;
    bool triangle_intersect(const Ray& ray, HitRecord& hit) const;

};

struct AABB
{
    glm::vec3 l;
    glm::vec3 u;

    AABB();
    AABB(const std::vector<Triangle>& triangles);
    AABB(const AABB b1, const AABB b2);

    uint8_t longest_axis() const;
    float surface_area() const;
    void update_box(const std::vector<Triangle>& triangles,int start,int end);
    void expand(Triangle t);
    bool box_intersect(const Ray& ray) const;
    void print() const
    {
    std::cout << "lower:" << l.x << "," << l.y << "," << l.z << std::endl;
    std::cout << "upper:" << u.x << "," << u.y << "," << u.z << std::endl;
    }
};



struct Rectangle {
    glm::vec3 origin;
    glm::vec3 w;
    glm::vec3 h;
};