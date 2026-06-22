#pragma once
#include "Geometry.hpp"
#include "tiny_obj_loader.h"
#include "Material.hpp"

struct Texture 
{
    int width, height, channels;
    float* data;

    Texture();
    Texture(std::string filename);
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
    bool load_object_from_file(std::string filename, std::string directory);
    Material make_default_material();
};