#include <SFML/Graphics.hpp>
#include <iostream>
#define GLM_ENABLE_EXPERIMENTAL
#define STB_IMAGE_IMPLEMENTATION
#include <tinycolormap.hpp>
#include <chrono>
#include "Geometry.hpp"
#include "Model.hpp"
#include "Scene.hpp"
#include "BVH.hpp"
#include "Renderer.hpp"
#include "Postprocessing.hpp"

using namespace glm;

constexpr bool HEATMAP = false;
#define LOG(code) \
    do { if constexpr (HEATMAP) { code; } } while(0)
struct HeatMap {
    uint8_t intersects;
    std::vector<uint8_t> heatpixels = {};
};
inline thread_local HeatMap heatMap;


std::vector<uint8_t> generateHeatmap(const std::vector<uint8_t>& pixels) {
    std::vector<uint8_t> heatmap;
    heatmap.reserve(pixels.size()*4);

    uint8_t minPixel = *std::min_element(pixels.begin(), pixels.end());
    uint8_t maxPixel = *std::max_element(pixels.begin(), pixels.end());
    float range = (maxPixel > minPixel) ? float(maxPixel - minPixel) : 1.0f;

    for (uint8_t val : pixels) {
        float normalized = (float(val) - float(minPixel)) / range;
        tinycolormap::Color c = tinycolormap::GetColor(normalized, tinycolormap::ColormapType::Turbo);

        heatmap.push_back(static_cast<uint8_t>(c.r() * 255.0f));
        heatmap.push_back(static_cast<uint8_t>(c.g() * 255.0f));
        heatmap.push_back(static_cast<uint8_t>(c.b() * 255.0f));
        heatmap.push_back(255);
    }
    return heatmap;
}


int main()
{ 
    /* ================= MODELLOADER ================= */

    std::cout << "[INFO] Loading models." << std::endl;

    std::vector<Model> models;
    std::vector<std::string> models_to_load = {"plane", "trees9"};
    for (const std::string& m : models_to_load)
        models.emplace_back(m, "/home/lukas/simple-raytracer/scene/");

    /* ================= SCENE SETUP ================= */

    std::cout << "[INFO] Initializing scene...";

    Scene scene(models, "/home/lukas/simple-raytracer/scene/setup.yaml");
    scene.update_data();

    std::cout << "Scene initialized (" << scene.triangles.size() << " triangles)." << std::endl;

    /* ================= BVH BUILDING ================= */

    std::cout << "[INFO] Building BVH...";
    auto bvh_start = std::chrono::high_resolution_clock::now();

    BVH bvh(scene.triangles);

    auto bvh_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> bvh_time = bvh_end - bvh_start;

    std::cout << "BVH construction finished in " << bvh_time.count() << " s." << std::endl;

    /* ================= RENDERING ================= */

    std::cout << "[INFO] Starting rendering...";
    auto render_start = std::chrono::high_resolution_clock::now();

    Renderer renderer(scene, bvh);
    auto pixel_values = renderer.render();

    auto render_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> render_time = render_end - render_start;
    std::cout << "Rendering completed in " << render_time.count() << " s." << std::endl;

    /* ================= POST-PROCESSING ================= */

    std::cout << "[INFO] Applying tone mapping...";
    auto pixel_colors = tone_mapping(pixel_values, 1.6);
    
    std::cout << "Pipeline finished successfully!" << std::endl;

    auto window = sf::RenderWindow(sf::VideoMode({scene.camera.resolution.x, scene.camera.resolution.y}), "RaytracerPOG");
    window.setFramerateLimit(60);
    

    sf::Texture texture(sf::Vector2u(scene.camera.resolution.x, scene.camera.resolution.y));
    LOG(texture.update(generateHeatmap(heatMap.heatpixels).data()););
    if (HEATMAP==false) texture.update(pixel_colors.data());

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
