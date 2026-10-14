#include <SFML/Graphics.hpp>
#include <iostream>
#include <filesystem>
#define GLM_ENABLE_EXPERIMENTAL
#define TINYGLTF_IMPLEMENTATION
#include <tiny_gltf.h>
#include <chrono>
#include "Geometry.hpp"
#include "Model.hpp"
#include "Scene.hpp"
#include "BVH.hpp"
#include "Renderer.hpp"
#include "Postprocessing.hpp"
#include "logger.hpp"

using namespace glm;
    
int main()
{ 
/* ================= SCENE SETUP ================= */

    std::cout << "[INFO] Initializing scene..." << std::endl;

    Scene scene("/home/lukas/simple-raytracer/scene/setup.yaml");
    scene.update_data();

    std::cout << "[INFO] Scene initialized (" << scene.triangles.size() << " triangles)." << std::endl;

    /* ================= BVH BUILDING ================= */

    std::cout << "[INFO] Building BVH..." << std::endl;
    auto bvh_start = std::chrono::high_resolution_clock::now();

    BVH bvh(scene.triangles);

    auto bvh_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> bvh_time = bvh_end - bvh_start;

    std::cout << "[INFO] BVH construction finished in " << bvh_time.count() << " s." << std::endl;

    /* ================= RENDERING ================= */

    std::cout << "[INFO] Starting rendering..." << std::endl;
    auto render_start = std::chrono::high_resolution_clock::now();

    Renderer renderer(scene, bvh);
    auto pixel_values = renderer.render();

    auto render_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> render_time = render_end - render_start;
    std::cout << "[INFO] Rendering completed in " << render_time.count() << " s." << std::endl;

    /* ================= POST-PROCESSING ================= */

    std::cout << "[INFO] Applying tone mapping..." << std::endl;
    auto pixel_colors = tone_mapping(pixel_values);
    
    std::cout << "[INFO] Pipeline finished successfully!" << std::endl;
    LOG(std::cout << "Total Intersection Tests:" << INT_COUNT/2100000 << std::endl;);

    auto window = sf::RenderWindow(sf::VideoMode({scene.camera.resolution.x, scene.camera.resolution.y}), "Raytracer");
    window.setFramerateLimit(60);

    std::filesystem::path dir = std::filesystem::path(PROJECT_ROOT) / "out/images";

    const std::string prefix = "output_";
    const std::string suffix = ".png";
    int max_index = 0;

    for (const auto& entry : std::filesystem::directory_iterator(dir))
    {
        auto name = entry.path().filename().string();

        if (name.starts_with(prefix) && name.ends_with(suffix))
        {
            auto number = name.substr(prefix.size(),name.size() - prefix.size() - suffix.size());
            max_index = std::max(max_index, std::stoi(number));
        }
    }

    sf::Image image(sf::Vector2u(scene.camera.resolution.x, scene.camera.resolution.y), pixel_colors.data());
    std::filesystem::path output = dir / ("output_" + std::to_string(max_index + 1) + ".png");
    if (!image.saveToFile(output))
    {
        throw std::runtime_error(
            "Failed to save image to " + output.string()
        );
    }

    sf::Texture texture(sf::Vector2u(scene.camera.resolution.x, scene.camera.resolution.y));
    LOG(texture.update(generateHeatmap(heat_map.heatpixels).data()););
    if (DEBUG==false) texture.update(pixel_colors.data());

    sf::Sprite sprite(texture);

    while (window.isOpen())
    {
        while (const std::optional event = window.pollEvent())
        {
            if (event->is<sf::Event::Closed>())
            {
                window.close();
            }
        }
        
        window.clear();
        window.draw(sprite);
        window.display();
    }
}
