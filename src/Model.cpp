#include "Model.hpp"
#include "stb_image.h"
#include <tiny_gltf.h>
#include <yaml-cpp/yaml.h>
#include <array>
#include <memory>
#include <cmath>


// Exact sRGB EOTF (IEC 61966-2-1), precomputed for all 8-bit values.
static float srgb8_to_linear(unsigned char c)
{
    static const std::array<float, 256> lut = [] {
        std::array<float, 256> t{};
        for (int i = 0; i < 256; ++i) {
            float s = i / 255.0f;
            t[i] = (s <= 0.04045f) ? s / 12.92f : std::pow((s + 0.055f) / 1.055f, 2.4f);
        }
        return t;
    }();
    return lut[c];
}

static int wrap_repeat(int i, int n)
{
    int r = i % n;
    return r < 0 ? r + n : r;
}


Texture::Texture() {}

Texture::Texture(std::string filename)
{
    int file_channels = 0;
    unsigned char* pixels = stbi_load(filename.c_str(), &width, &height, &file_channels, 3);
    if (!pixels) {
        std::cerr << "ERR: could not load texture " << filename << std::endl;
        width = height = 0;
        return;
    }
    load_srgb8(width, height, 3, pixels);
    stbi_image_free(pixels);
}

Texture::Texture(int w, int h, int c, const unsigned char* in_data)
{
    load_srgb8(w, h, c, in_data);
}

void Texture::load_srgb8(int w, int h, int c, const unsigned char* in_data)
{
    width = w;
    height = h;
    std::shared_ptr<float[]> pixels(new float[width * height * channels]);
    
    for (int i = 0; i < width * height; ++i) {
        const unsigned char* px = in_data + i * c;
        unsigned char r = px[0];
        unsigned char g = (c >= 3) ? px[1] : px[0];
        unsigned char b = (c >= 3) ? px[2] : px[0];

        pixels[i * channels + 0] = srgb8_to_linear(r);
        pixels[i * channels + 1] = srgb8_to_linear(g);
        pixels[i * channels + 2] = srgb8_to_linear(b);
    }
    data = std::move(pixels);
}

glm::vec3 Texture::sample(glm::vec2 uv) const
{
    if (!data || width <= 0 || height <= 0)
        return glm::vec3(1.0f);

    float x = uv.x * width  - 0.5f;
    float y = uv.y * height - 0.5f;
    float fx = std::floor(x);
    float fy = std::floor(y);
    float tx = x - fx;
    float ty = y - fy;

    int x0 = wrap_repeat(int(fx), width);
    int y0 = wrap_repeat(int(fy), height);
    int x1 = wrap_repeat(x0 + 1, width);
    int y1 = wrap_repeat(y0 + 1, height);

    auto texel = [&](int xi, int yi) {
        const float* p = data.get() + (yi * width + xi) * channels;
        return glm::vec3(p[0], p[1], p[2]);
    };

    glm::vec3 top    = glm::mix(texel(x0, y0), texel(x1, y0), tx);
    glm::vec3 bottom = glm::mix(texel(x0, y1), texel(x1, y1), tx);
    return glm::mix(top, bottom, ty);
}

Model::Model(std::string filename, std::string directory)
{
    load_object_from_file(filename, directory);
}

int Model::load_texture(const tinygltf::Image& image) 
{
    if (image.image.empty()) return -1;
    textures.emplace_back(image.width, image.height, image.component, image.image.data());
    return static_cast<int>(textures.size() - 1);
}

bool Model::load_object_from_file(std::string filename, std::string directory)
{
    name = filename;
    tinygltf::Model gltf_model;
    tinygltf::TinyGLTF loader;
    std::string err;
    std::string warn;
    
    std::string filepath = directory + filename + "/" + filename + ".glb";
    bool ret = loader.LoadBinaryFromFile(&gltf_model, &err, &warn, filepath);
    
    if (!warn.empty()) std::cout << "WARN: " << warn << std::endl;
    if (!err.empty()) {
        std::cerr << "ERR: " << err << std::endl;
        return false;
    }
    if (!ret) return false;

    materials.reserve(gltf_model.materials.size() + 1);
    materials.push_back(make_default_material());

    for (const auto& mat : gltf_model.materials) {
        Material m;
        if (mat.pbrMetallicRoughness.baseColorFactor.size() >= 3) 
            m.albedo = glm::vec3(mat.pbrMetallicRoughness.baseColorFactor[0], mat.pbrMetallicRoughness.baseColorFactor[1], mat.pbrMetallicRoughness.baseColorFactor[2]);
        
        m.metallic = mat.pbrMetallicRoughness.metallicFactor;
        m.roughness = mat.pbrMetallicRoughness.roughnessFactor;
        
        if (mat.emissiveFactor.size() >= 3)
            m.emissive = 10.f * glm::vec3(mat.emissiveFactor[0], mat.emissiveFactor[1], mat.emissiveFactor[2]);
        
        auto ext_it = mat.extensions.find("KHR_materials_ior");
        if (ext_it != mat.extensions.end()) {
            if (ext_it->second.Has("ior")) {
                m.ior = ext_it->second.Get("ior").GetNumberAsDouble();
            }
        } else {
            m.ior = 1.45f;
        }

        int tex_id = mat.pbrMetallicRoughness.baseColorTexture.index;
        if (tex_id >= 0 && tex_id < gltf_model.textures.size()) {
            int image_id = gltf_model.textures[tex_id].source;
            if (image_id >= 0) m.diffuseTex = load_texture(gltf_model.images[image_id]);
        }
        materials.push_back(m);
    }

    size_t total_faces = 0;
    for (const auto& mesh : gltf_model.meshes) {
        for (const auto& primitive : mesh.primitives) {
            if (primitive.mode != TINYGLTF_MODE_TRIANGLES) continue;
            if (primitive.indices >= 0) {
                total_faces += gltf_model.accessors[primitive.indices].count / 3;
            } else {
                if (primitive.attributes.find("POSITION") != primitive.attributes.end()) {
                    total_faces += gltf_model.accessors[primitive.attributes.at("POSITION")].count / 3;
                }
            }
        }
    }
    triangles.reserve(total_faces);

    for (const auto& mesh : gltf_model.meshes) {
        for (const auto& primitive : mesh.primitives) {
            if (primitive.mode != TINYGLTF_MODE_TRIANGLES) {
                std::cout << "Warning: non-triangle primitive skipped." << std::endl;
                continue;
            }

            const float* positions = nullptr;
            const float* normals = nullptr;
            const float* uvs = nullptr;
            
            int posStride = 0, normStride = 0, uvStride = 0;
            
            if (primitive.attributes.find("POSITION") != primitive.attributes.end()) {
                const tinygltf::Accessor& accessor = gltf_model.accessors[primitive.attributes.at("POSITION")];
                const tinygltf::BufferView& view = gltf_model.bufferViews[accessor.bufferView];
                positions = reinterpret_cast<const float*>(&gltf_model.buffers[view.buffer].data[view.byteOffset + accessor.byteOffset]);
                posStride = accessor.ByteStride(view);
            }
            if (primitive.attributes.find("NORMAL") != primitive.attributes.end()) {
                const tinygltf::Accessor& accessor = gltf_model.accessors[primitive.attributes.at("NORMAL")];
                const tinygltf::BufferView& view = gltf_model.bufferViews[accessor.bufferView];
                normals = reinterpret_cast<const float*>(&gltf_model.buffers[view.buffer].data[view.byteOffset + accessor.byteOffset]);
                normStride = accessor.ByteStride(view);
            }
            if (primitive.attributes.find("TEXCOORD_0") != primitive.attributes.end()) {
                const tinygltf::Accessor& accessor = gltf_model.accessors[primitive.attributes.at("TEXCOORD_0")];
                const tinygltf::BufferView& view = gltf_model.bufferViews[accessor.bufferView];
                uvs = reinterpret_cast<const float*>(&gltf_model.buffers[view.buffer].data[view.byteOffset + accessor.byteOffset]);
                uvStride = accessor.ByteStride(view);
            }

            int mat_id = primitive.material >= 0 ? primitive.material + 1 : 0;
            
            if (primitive.indices >= 0) {
                const tinygltf::Accessor& indexAccessor = gltf_model.accessors[primitive.indices];
                const tinygltf::BufferView& indexBufferView = gltf_model.bufferViews[indexAccessor.bufferView];
                const tinygltf::Buffer& indexBuffer = gltf_model.buffers[indexBufferView.buffer];
                const unsigned char* indexData = &indexBuffer.data[indexBufferView.byteOffset + indexAccessor.byteOffset];
                int indexStride = indexAccessor.ByteStride(indexBufferView);
                
                for (size_t i = 0; i < indexAccessor.count; i += 3) {
                    Triangle tri;
                    tri.mat_id = mat_id;
                    tri.hasVertexNormals = (normals != nullptr);
                    
                    for (int j = 0; j < 3; ++j) {
                        uint32_t idx = 0;
                        if (indexAccessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT) {
                            idx = *reinterpret_cast<const uint16_t*>(indexData + (i + j) * indexStride);
                        } else if (indexAccessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT) {
                            idx = *reinterpret_cast<const uint32_t*>(indexData + (i + j) * indexStride);
                        } else {
                            idx = *(indexData + (i + j) * indexStride);
                        }
                        
                        if (positions) {
                            const float* p = reinterpret_cast<const float*>(reinterpret_cast<const unsigned char*>(positions) + idx * posStride);
                            tri.v[j] = glm::vec3(p[0], p[1], p[2]);
                        }
                        if (normals) {
                            const float* n = reinterpret_cast<const float*>(reinterpret_cast<const unsigned char*>(normals) + idx * normStride);
                            tri.n[j] = glm::vec3(n[0], n[1], n[2]);
                        }
                        if (uvs) {
                            const float* u = reinterpret_cast<const float*>(reinterpret_cast<const unsigned char*>(uvs) + idx * uvStride);
                            tri.uv[j] = glm::vec2(u[0], u[1]);
                        }
                    }
                    triangles.push_back(tri);
                }
            } else {
                if (!positions) continue;
                const tinygltf::Accessor& accessor = gltf_model.accessors[primitive.attributes.at("POSITION")];
                for (size_t i = 0; i < accessor.count; i += 3) {
                    Triangle tri;
                    tri.mat_id = mat_id;
                    tri.hasVertexNormals = (normals != nullptr);
                    for (int j = 0; j < 3; ++j) {
                        uint32_t idx = i + j;
                        if (positions) {
                            const float* p = reinterpret_cast<const float*>(reinterpret_cast<const unsigned char*>(positions) + idx * posStride);
                            tri.v[j] = glm::vec3(p[0], p[1], p[2]);
                        }
                        if (normals) {
                            const float* n = reinterpret_cast<const float*>(reinterpret_cast<const unsigned char*>(normals) + idx * normStride);
                            tri.n[j] = glm::vec3(n[0], n[1], n[2]);
                        }
                        if (uvs) {
                            const float* u = reinterpret_cast<const float*>(reinterpret_cast<const unsigned char*>(uvs) + idx * uvStride);
                            tri.uv[j] = glm::vec2(u[0], u[1]);
                        }
                    }
                    triangles.push_back(tri);
                }
            }
        }
    }
    return true;
}

bool Model::load_object_from_file_legacy(std::string filename, std::string directory)
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
        mat.emissive = {m.emission[0], m.emission[1], m.emission[2]};

        glm::vec3 ks(m.specular[0], m.specular[1], m.specular[2]);
        glm::vec3 kd(m.diffuse[0], m.diffuse[1], m.diffuse[2]);

        float ks_lum = 0.2126f * ks.r + 0.7152f * ks.g + 0.0722f * ks.b;
        float kd_lum = 0.2126f * kd.r + 0.7152f * kd.g + 0.0722f * kd.b;

        if (ks_lum > 0.01f) {
            mat.metallic = 1.0f;
            mat.albedo = ks;
        } else {
            mat.metallic = 0.0f;
            mat.albedo = kd;
        }
        if (m.shininess > 0.0f) {
            mat.roughness = glm::sqrt(2.0f / (m.shininess + 2.0f));
        } else {
            mat.roughness = 1.0f;
        }
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