#pragma once
#include <tinycolormap.hpp>

inline constexpr bool DEBUG = true;
inline int INT_COUNT = 0;

#define LOG(code) do { if constexpr (DEBUG) { code; } } while(0)


struct HeatMap {
    uint8_t intersects = 0;
    std::vector<uint8_t> heatpixels = {};
};

inline HeatMap heat_map;

inline std::vector<uint8_t> generateHeatmap(const std::vector<uint8_t>& pixels) {
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