#pragma once
#include <vector>
#include <string>
#include <unordered_set>
#include "Geometry.hpp"
#include "Model.hpp"

struct Object {
    const Model* model;
    glm::mat4 transform;

    Object(const std::string& model_name, const glm::mat4& transform, const std::vector<Model>& models);
};


struct LightSource {
    Rectangle area;
    glm::vec3 color;
    int count;

    LightSource(const glm::vec3& origin, const glm::vec3& w, const glm::vec3& h, const glm::vec3& color, int count);

    std::vector<glm::vec3> get_lightpoints() const;
};


struct Camera {
    glm::vec3 position;
    glm::vec3 forward;
    glm::vec3 up;
    float fov;
    glm::uvec2 resolution;
    Rectangle image_plane;
    
    Camera();
    Camera(const glm::vec3& position, const glm::vec3& forward, const glm::vec3& up, float fov, int width_pixels, int height_pixels);

    std::vector<Ray> generate_rays() const;
    Ray generate_ray(int px, int py) const;


private: 
    void build_image_plane(float distance_to_camera);
};


class Scene {
public:
    std::vector<Model> models;
    std::vector<Object> objects;

    Camera camera;
    std::vector<LightSource> lights;

    std::vector<Triangle> triangles;
    std::vector<Material> materials;
    std::vector<Texture> textures;

    Scene(const std::string& filename);
    void update_data();

private:
    void setup_scene(const std::string& filename);
    std::unordered_set<std::string> get_model_names(const std::string& filename);
};