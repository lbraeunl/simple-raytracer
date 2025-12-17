#pragma once
#include <vector>
#include <string>
#include "Geometry.hpp"
#include "Model.hpp"

struct Object {
    const Model* model;
    glm::mat4 transform;

    Object(const std::string& model_name, const glm::mat4& transform, const std::vector<Model>& models);
};


struct LightSource {
    glm::vec3 position;
    glm::vec3 color;

    LightSource(const glm::vec3& position, const glm::vec3& color);
};


struct Camera {
    glm::vec3 position;
    glm::vec3 forward;
    glm::vec3 up;
    float fov;
    glm::ivec2 resolution;
    Rectangle image_plane;
    
    Camera();
    Camera(const glm::vec3& position, const glm::vec3& forward, const glm::vec3& up, float fov, int width_pixels, int height_pixels);

private: 
    void build_image_plane(float distance_to_camera);
};


class Scene {
public:
    Camera camera;
    std::vector<Object> objects;
    std::vector<LightSource> lights;

    std::vector<Triangle> triangles;
    std::vector<Material> materials;
    std::vector<Texture> textures;

    Scene(const std::vector<Model>& models, const std::string& filename);
    void update_data();

private:
    void setup_scene(const std::vector<Model>& models, const std::string& filename);

};