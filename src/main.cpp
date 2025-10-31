#include <SFML/Graphics.hpp>
#include "tiny_obj_loader.h"
#include <iostream>
#include <glm/vec3.hpp>

using namespace glm;

// bool load_from_obj(const char* filename)
// {
    
// }

bool is_there_something(vec3 direction, auto shapes, auto attrib)
{
    for (const auto& shape : shapes) {
        size_t index_offset = 0;

        for (size_t f = 0; f < shape.mesh.num_face_vertices.size(); f++) {
            int fv = shape.mesh.num_face_vertices[f];  // Anzahl Vertices pro Face (meist 3)

            for (size_t v = 0; v < fv; v++) {
                // Hole den Index des Vertex
                tinyobj::index_t idx = shape.mesh.indices[index_offset + v];

                // Zugriff auf die eigentlichen Vertexkoordinaten
                float vx = attrib.vertices[3 * idx.vertex_index + 0];
                float vy = attrib.vertices[3 * idx.vertex_index + 1];
                float vz = attrib.vertices[3 * idx.vertex_index + 2];

                std::cout << "Face " << f << ", Vertex " << v << ": "
                        << vx << ", " << vy << ", " << vz << std::endl;
            }

            index_offset += fv;  // Weiter zum nächsten Face
        }
    }
}

int main()
{

    vec3 cam = {0.f,-40.f,10.f};
    vec3 img_point = {-5.f,-30.f,5.f};
    vec3 veci = {10.f,0.f,0.f};
    vec3 vecj = {0.f,0.f,10.f};
    int width = 128;
    int height = 128;
    std::vector<std::vector<bool>> grid(height, std::vector<bool>(width, false));
    

    //attrib will contain the vertex arrays of the file
	tinyobj::attrib_t attrib;
    //shapes contains the info for each separate object in the file
	std::vector<tinyobj::shape_t> shapes;
    //materials contains the information about the material of each shape, but we won't use it.
    std::vector<tinyobj::material_t> materials;

    //error and warning output from the load function
	std::string warn;
	std::string err;

    //load the OBJ file
	tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, "./tinker.obj", nullptr);

    // for (size_t i = 0; i < attrib.vertices.size(); i += 3) {
    //     float x = attrib.vertices[i + 0];
    //     float y = attrib.vertices[i + 1];
    //     float z = attrib.vertices[i + 2];

    //     std::cout << "Vertex " << (i / 3) << ": "
    //             << x << ", " << y << ", " << z << std::endl;
    // }

    



    //make sure to output the warnings to the console, in case there are issues with the file
	if (!warn.empty()) {
		std::cout << "WARN: " << warn << std::endl;
	}
    //if we have any error, print it to the console, and break the mesh loading.
    //This happens if the file can't be found or is malformed
	if (!err.empty()) {
		std::cerr << err << std::endl;
		return false;
	}

    for(int i=0;i<width;++i) 
    {
        for(int j=0;j<height;++j) 
        {
            vec3 pos = img_point+(veci*float(i))/float(width)+(vecj*float(j))/float(height);
            vec3 dir = pos-cam;
        }
    }

    // auto window = sf::RenderWindow(sf::VideoMode({1920u, 1080u}), "CMake SFML Project");
    // window.setFramerateLimit(144);

    // while (window.isOpen())
    // {
    //     while (const std::optional event = window.pollEvent())
    //     {
    //         if (event->is<sf::Event::Closed>())
    //         {
    //             window.close();
    //         }
    //     }

    //     window.clear();
    //     window.display();
    // }
}
