#include "Scene.hpp"
#include <yaml-cpp/yaml.h>
#include <glm/gtc/matrix_transform.hpp>


Object::Object(const std::string& model_name,const glm::mat4& transform, const auto& models):transform(transform)
{
    model = &models.at(model_name);
}


LightSource::LightSource(const glm::vec3& position, const glm::vec3& color): position(position), color(color){}


Camera::Camera(const glm::vec3& position, const glm::vec3& forward, const glm::vec3& up, float fov, int width_pixels, int height_pixels): position(position), forward(forward), up(up), fov(fov)
{
    resolution.x = width_pixels;
    resolution.y = height_pixels;
    build_image_plane(0.01f);
}


void Camera::build_image_plane(float distance_to_camera)
{
    glm::vec3 right = glm::normalize(glm::cross(forward, up));
    up = glm::normalize(glm::cross(right, forward));

    float height = 2.0f * distance_to_camera * tan(glm::radians(fov) * 0.5f);
    float aspect = float(resolution.x) / float(resolution.y);
    float width = height * aspect;

    image_plane.h = up*height;
    image_plane.w = right*width;
    image_plane.origin = position + forward*distance_to_camera - 0.5f*image_plane.h - 0.5f*image_plane.w;
}

Scene::Scene(const auto& models, const std::string& filename)
{
    setup_scene(models, filename);
}


void Scene::setup_scene( const auto& models, const std::string& filename)
{
    YAML::Node file = YAML::LoadFile(filename);

    //Camera
    const auto& camNode = file["camera"];
    camera = Camera(
        glm::vec3(camNode["position"][0].as<float>(), camNode["position"][1].as<float>(), camNode["position"][2].as<float>()),
        glm::vec3(camNode["forward"][0].as<float>(), camNode["forward"][1].as<float>(), camNode["forward"][2].as<float>()),
        glm::vec3(camNode["up"][0].as<float>(), camNode["up"][1].as<float>(), camNode["up"][2].as<float>()),
        camNode["fov"].as<float>(),
        camNode["resolution"][0].as<int>(),
        camNode["resolution"][1].as<int>()
    );

    //Lights
    for (const auto& lightNode : file["lights"]) 
    {
        glm::vec3 pos(lightNode["position"][0].as<float>(),lightNode["position"][1].as<float>(),lightNode["position"][2].as<float>());
        glm::vec3 color(lightNode["color"][0].as<float>(),lightNode["color"][1].as<float>(),lightNode["color"][2].as<float>());
        lights.emplace_back(pos, color);
    }

    //Objects
    for (const auto& objNode : file["objects"]) 
    {
        glm::mat4 transform(1.0f);
        if (objNode["scale"]) {
            auto s = objNode["scale"];
            glm::vec3 scale(s[0].as<float>(), s[1].as<float>(), s[2].as<float>());
            transform = glm::scale(transform, scale);
        }
        if (objNode["rotate"]) {
            auto r = objNode["rotate"];
            glm::vec3 rotDeg(r[0].as<float>(), r[1].as<float>(), r[2].as<float>());
            glm::vec3 rot = glm::radians(rotDeg);

            transform = glm::rotate(transform, rot.x, glm::vec3(1,0,0));
            transform = glm::rotate(transform, rot.y, glm::vec3(0,1,0));
            transform = glm::rotate(transform, rot.z, glm::vec3(0,0,1));
        }
        if (objNode["translate"]) {
            auto t = objNode["translate"];
            glm::vec3 pos(t[0].as<float>(), t[1].as<float>(), t[2].as<float>());
            transform = glm::translate(transform, pos);
        }
        objects.emplace_back(file["model_name"].as<std::string>(),transform, models);
    }

}