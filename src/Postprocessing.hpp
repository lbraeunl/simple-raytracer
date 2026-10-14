#pragma once
#include "Geometry.hpp"

// HDR linear radiance -> 8-bit sRGB RGBA
std::vector<glm::uint8> tone_mapping(const std::vector<glm::vec3>& pixels_in);
