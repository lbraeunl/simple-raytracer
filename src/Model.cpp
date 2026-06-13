#include "Model.hpp"
#include "stb_image.h"
#include <yaml-cpp/yaml.h>



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

    materials.reserve(inMaterials.size()+1);
    materials.push_back(make_default_material());

    for (auto& m : inMaterials) {
        Material mat;

        mat.albedo = {m.diffuse[0], m.diffuse[1], m.diffuse[2]};
        mat.emissive = {m.emission[0], m.emission[1], m.emission[2]};
        mat.roughness = m.shininess;
        // if (m.shininess > 0.0f) {
        //     mat.roughness = glm::sqrt(2.0f / (m.shininess + 2.0f));
        // } else {
        //     mat.roughness = 1.0f;
        // }
        mat.metallic = (m.specular[0] > 0.0f) ? 1.0f : 0.0f; 
        mat.ior = m.ior;

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

    for (const auto& shape : inShapes) 
    {
        size_t index_offset = 0;
        for (size_t f = 0; f < shape.mesh.num_face_vertices.size(); f++) 
        {
            size_t fv = shape.mesh.num_face_vertices[f];
            if (fv != 3) 
            {
                std::cout << "Warning:" << filename << " contains a non-triangular shape that is skipped." << std::endl;
                index_offset += fv;
                continue;
            }

            Triangle tria;
            for(int i=0;i<3;++i)
            {
                tinyobj::index_t idx = shape.mesh.indices[index_offset+i];
                tria.v[i] = glm::vec3(inAttrib.vertices[3*idx.vertex_index+0], inAttrib.vertices[3*idx.vertex_index+1], inAttrib.vertices[3*idx.vertex_index+2]);

                if (idx.normal_index >= 0)
                {
                    tria.n[i] = glm::vec3(inAttrib.normals[3*idx.normal_index + 0], inAttrib.normals[3*idx.normal_index + 1], inAttrib.normals[3*idx.normal_index + 2]);
                    tria.hasVertexNormals = true;
                }

                if (idx.texcoord_index >= 0)
                {
                    tria.uv[i]= glm::vec2(inAttrib.texcoords[2*idx.texcoord_index+0],inAttrib.texcoords[2*idx.texcoord_index+1]);
                }
            }
            
            tria.mat_id = shape.mesh.material_ids[f]+1;
            triangles.push_back(tria);

            index_offset += 3;
        }
    }
    return true;
}

int Model::load_texture(std::string tex_directory, std::string texname) 
{
    if (texname.empty()) return -1;

    textures.emplace_back(tex_directory + "/" + texname);

    return static_cast<int>(textures.size() - 1);
}

Material Model::make_default_material()
{
    Material mat;

    mat.albedo = glm::vec3(0.7f);
    mat.emissive = glm::vec3(0.0f);
    mat.metallic = 0.0f;
    mat.roughness = 0.5f;
    mat.ior = 1.45f;

    return mat;
}