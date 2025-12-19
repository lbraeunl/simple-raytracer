#pragma once
#include "Geometry.hpp"

inline float reinhard(float x);

std::vector<glm::uint8> tone_mapping(const std::vector<glm::vec3>& pixel_floats, float yamma);