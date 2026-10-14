#pragma once
#include "Geometry.hpp"
#include "tiny_obj_loader.h"
#include "Material.hpp"
#include <memory>

namespace tinygltf {
    struct Image;
}


struct Texture
{
    int width = 0, height = 0;
    static constexpr int channels = 3;
    
    std::shared_ptr<const float[]> data;

    Texture();
    Texture(std::string filename);
    Texture(int w, int h, int c, const unsigned char* in_data);

    glm::vec3 sample(glm::vec2 uv) const;

private:
    void load_srgb8(int w, int h, int c, const unsigned char* in_data);
};



class Model {
public:
    std::vector<Material> materials;
    std::vector<Texture> textures;
    std::vector<Triangle> triangles;
    std::string name;

    Model(std::string filename, std::string directory);
    
private:
    int load_texture(std::string tex_directory, std::string texname);
    int load_texture(const tinygltf::Image& image);
    bool load_object_from_file(std::string filename, std::string directory);
    bool load_object_from_file_legacy(std::string filename, std::string directory);
    Material make_default_material();
};