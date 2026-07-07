#pragma once
#include "Geometry.hpp"
#include "tiny_obj_loader.h"
#include "Material.hpp"

namespace tinygltf {
    struct Image;
}

struct Texture 
{
    int width, height, channels;
    float* data;

    Texture();
    Texture(std::string filename);
    Texture(int w, int h, int c, const unsigned char* in_data);
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