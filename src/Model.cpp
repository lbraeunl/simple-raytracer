#include "Model.hpp"
#include <unordered_map>
#include "stb_image.h"


Texture::Texture(std::string filename)
{
    data = stbi_loadf(filename.c_str(), &width, &height, &channels, 3);
}

Model::Model(std::string filename, std::string directory)
{
    load_object_from_file(filename, directory);
}

bool Model::load_object_from_file(std::string filename, std::string directory)
{
    name = filename;
    tinyobj::attrib_t inAttrib;
    std::vector<tinyobj::shape_t> inShapes;
    std::vector<tinyobj::material_t> inMaterials;
    std::string warn, err;

    bool ret = tinyobj::LoadObj(&inAttrib, &inShapes, &inMaterials,&warn, &err,(directory+filename+"/"+filename+".obj").c_str(),(directory+filename+"/").c_str());

    if (!warn.empty()) {
		std::cout << "WARN: " << warn << std::endl;
	}

	if (!err.empty()) {
		std::cerr << err << std::endl;
		return false;
	}

    materials.reserve(inMaterials.size());
    for (auto& m : inMaterials) {
        Material mat;

        mat.ambientColor = {m.ambient[0], m.ambient[1], m.ambient[2]};
        mat.diffuseColor = {m.diffuse[0], m.diffuse[1], m.diffuse[2]};
        mat.specularColor = {m.specular[0], m.specular[1], m.specular[2]};

        mat.ambientTex  = load_texture(directory+"/" + filename, m.ambient_texname);
        mat.diffuseTex  = load_texture(directory+"/" + filename, m.diffuse_texname);
        mat.specularTex = load_texture(directory+"/" + filename, m.specular_texname);

        materials.push_back(mat);
    }

    size_t total_faces = 0;
    for (const auto& shape : inShapes) {
        total_faces += shape.mesh.num_face_vertices.size();
    }

    triangles.reserve(total_faces);

    Triangle tria;
    for (const auto& shape : inShapes) 
    {
        for (size_t f = 0; f < shape.mesh.num_face_vertices.size(); f++) 
        {
            if (shape.mesh.num_face_vertices[f] != 3) {
                std::cout << "Warning:" << filename << " contains a non-triangular shape that is skipped." << std::endl;
                continue;
            }

            for(int i=0;i<3;++i)
            {
                tinyobj::index_t idx = shape.mesh.indices[3*f+i];
                tria.v[i] = glm::vec3(inAttrib.vertices[3*idx.vertex_index+0], inAttrib.vertices[3*idx.vertex_index+1], inAttrib.vertices[3*idx.vertex_index+2]);
                if (idx.texcoord_index >= 0)
                {
                    tria.uv[i]= glm::vec2(inAttrib.texcoords[2*idx.texcoord_index+0],inAttrib.texcoords[2*idx.texcoord_index+1]);
                }
            }
            tria.mat_id = shape.mesh.material_ids[f];
            triangles.emplace_back(tria);
        }
    }
    return true;
}

int Model::load_texture(std::string tex_directory, std::string texname) 
{
    if (texname.empty()) return -1;

    Texture tex(tex_directory + "/" + texname);
    textures.push_back(tex);

    return static_cast<int>(textures.size() - 1);
}

// void Model::add_floor(glm::vec3 color,int size)
// {
//     Material mat;
//     mat.diffuseColor = color;
//     materials.push_back(mat);
//     triangles.emplace_back(glm::vec3(size*0.5f,size*(-0.5f),0.f),glm::vec3(size*0.5f, size*0.5f,0.f),glm::vec3(size*-0.5f, size*-0.5f,0.f),int(materials.size()-1));
//     triangles.emplace_back(glm::vec3(size*0.5f, size*0.5f,0.f),glm::vec3(size*(-0.5f),size*0.5f,0.f),glm::vec3(size*-0.5f, size*-0.5f,0.f),int(materials.size()-1));
// }