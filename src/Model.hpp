#pragma once
#include <iostream>
#include "Geometry.hpp"
#include "tiny_obj_loader.h"

struct Texture 
{
    int width, height, channels;
    float* data;

    Texture();
    Texture(std::string filename);
};


struct Material 
{
    glm::vec3 ambientColor{1.0f};
    glm::vec3 diffuseColor{1.0f};
    glm::vec3 specularColor{1.0f};

    int ambientTex = -1;
    int diffuseTex = -1;
    int specularTex = -1;
};


class Model {
public:
    std::vector<Material> materials;
    std::vector<Texture> textures;
    std::vector<Triangle> triangles;
    std::string name;

    Model(std::string filename, std::string directory);
    //void add_floor(glm::vec3 color,int size);

private:
    int load_texture(std::string tex_directory, std::string texname);
    bool load_object_from_file(std::string filename, std::string directory);

};