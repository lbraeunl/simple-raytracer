#include "Scene.hpp"
#include <yaml-cpp/yaml.h>
#include <glm/gtc/matrix_transform.hpp>


Object::Object(const std::string& model_name,const glm::mat4& transform, const std::vector<Model>& models):transform(transform)
{
    for(const Model& m : models)
    {
        if (m.name == model_name) model = &m;
    }
}


LightSource::LightSource(const glm::vec3& position, const glm::vec3& color): position(position), color(color){}


Camera::Camera() : position(0.0f, 0.0f, 0.0f), forward(1.0f, 0.0f, 0.0f), up(0.0f, 1.0f, 0.0f), fov(60.0f), resolution(glm::vec2({1080,720})) {}

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

std::vector<Ray> Camera::generate_rays() const
{
    size_t pixel_count = resolution.x*resolution.y;
    std::vector<Ray> rays;
    rays.reserve(pixel_count);
    
    for (size_t i = 0; i < pixel_count; ++i)
    {
        float u = (i % resolution.x + 0.5f) / float(resolution.x);
        float v = (i / resolution.x + 0.5f) / float(resolution.y);

        glm::vec3 image_pos = image_plane.origin + u*image_plane.w + v*image_plane.h;
        rays.emplace_back(position, image_pos);
    }
    return rays;
}


Scene::Scene(const std::vector<Model>& models, const std::string& filename)
{
    setup_scene(models, filename);
}


void Scene::setup_scene(const std::vector<Model>& models, const std::string& filename)
{
    YAML::Node file = YAML::LoadFile(filename);

    const auto& camNode = file["camera"];
    camera = Camera(
        glm::vec3(camNode["position"][0].as<float>(), camNode["position"][1].as<float>(), camNode["position"][2].as<float>()),
        glm::vec3(camNode["forward"][0].as<float>(), camNode["forward"][1].as<float>(), camNode["forward"][2].as<float>()),
        glm::vec3(camNode["up"][0].as<float>(), camNode["up"][1].as<float>(), camNode["up"][2].as<float>()),
        camNode["fov"].as<float>(),
        camNode["resolution"][0].as<int>(),
        camNode["resolution"][1].as<int>()
    );

    for (const auto& lightNode : file["lights"]) 
    {
        glm::vec3 pos(lightNode["position"][0].as<float>(),lightNode["position"][1].as<float>(),lightNode["position"][2].as<float>());
        glm::vec3 color(lightNode["color"][0].as<float>(),lightNode["color"][1].as<float>(),lightNode["color"][2].as<float>());
        lights.emplace_back(pos, color);
    }

    for (const auto& objNode : file["objects"]) 
    {
        glm::mat4 transform(1.0f);

        if (objNode["scale"]) 
        {
            auto s = objNode["scale"];
            glm::vec3 scale(s[0].as<float>(), s[1].as<float>(), s[2].as<float>());
            transform = glm::scale(transform, scale);
        }

        if (objNode["rotate"]) 
        {
            auto r = objNode["rotate"];
            glm::vec3 rotDeg(r[0].as<float>(), r[1].as<float>(), r[2].as<float>());
            glm::vec3 rot = glm::radians(rotDeg);

            transform = glm::rotate(transform, rot.x, glm::vec3(1,0,0));
            transform = glm::rotate(transform, rot.y, glm::vec3(0,1,0));
            transform = glm::rotate(transform, rot.z, glm::vec3(0,0,1));
        }

        if (objNode["translate"]) 
        {
            auto t = objNode["translate"];
            glm::vec3 pos(t[0].as<float>(), t[1].as<float>(), t[2].as<float>());
            transform = glm::translate(transform, pos);
        }

        objects.emplace_back(objNode["model_name"].as<std::string>(),transform, models);
    }
    //TODO: Andere Objekte wie Boden, etc. implementieren
}


void Scene::update_data()
{
    size_t total_triangles = 0;
    size_t total_materials = 0;
    size_t total_textures = 0;

    for (const auto& obj : objects) {
        total_triangles += obj.model->triangles.size();
        total_materials += obj.model->materials.size();
        total_textures += obj.model->textures.size();
    }

    triangles.reserve(total_triangles);
    materials.reserve(total_materials);
    textures.reserve(total_textures);

    for (const auto& obj : objects)
    {      
        int material_offset = materials.size();
        int texture_offset  = textures.size();

        textures.insert(textures.end(), obj.model->textures.begin(), obj.model->textures.end());

        for (const auto& mat : obj.model->materials) 
        {
            Material m = mat;
            if (m.diffuseTex >= 0)
                m.diffuseTex += texture_offset;
            materials.push_back(m);
        }

        for (const auto& tri : obj.model->triangles) 
        {
            Triangle t;
            t.mat_id = tri.mat_id + material_offset;
            for (int i = 0; i < 3; i++) 
            {
                t.v[i] = glm::vec3(obj.transform * glm::vec4(tri.v[i], 1.0f));
                t.uv[i] = tri.uv[i];
            }
            t.update();
            triangles.push_back(t);
        }
    }
}