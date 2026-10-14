#include "Postprocessing.hpp"
#include <algorithm>
#include <cmath>


static float luminance(const glm::vec3& c)
{
    return 0.2126f * c.r + 0.7152f * c.g + 0.0722f * c.b;
}


// Reinhard et al. 2002, applied to luminance so hue and saturation are kept.
static glm::vec3 reinhard(const glm::vec3& c)
{
    float l_in = luminance(c);
    if (l_in <= 0.0f) return glm::vec3(0.0f);
    float l_out = l_in / (1.0f + l_in);
    return c * (l_out / l_in);
}


// Exact sRGB OETF (IEC 61966-2-1).
static float linear_to_srgb(float x)
{
    return (x <= 0.0031308f) ? 12.92f * x : 1.055f * std::pow(x, 1.0f / 2.4f) - 0.055f;
}


std::vector<glm::uint8> tone_mapping(const std::vector<glm::vec3>& pixels_in)
{
    std::vector<glm::uint8> pixels_out(pixels_in.size()*4);
    for(size_t i=0; i<pixels_in.size(); ++i)
    {
        glm::vec3 c = pixels_in[i];
        if (!std::isfinite(c.r) || !std::isfinite(c.g) || !std::isfinite(c.b))
            c = glm::vec3(0.0f);

        c = reinhard(glm::max(c, glm::vec3(0.0f)));

        for(int ch=0; ch<3; ++ch){
            float v = std::clamp(linear_to_srgb(c[ch]), 0.0f, 1.0f);
            pixels_out[4*i+ch] = glm::uint8(v * 255.0f + 0.5f);
        }
        pixels_out[4*i+3] = 255;
    }
    return pixels_out;
}
