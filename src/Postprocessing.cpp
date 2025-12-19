#include "Postprocessing.hpp"


inline float reinhard(float x)
{
    return x / (1.0f + x);
}


std::vector<glm::uint8> tone_mapping(const std::vector<glm::vec3>& pixels_in, float yamma) //extend to different functions
{
    std::vector<glm::uint8> pixels_out(pixels_in.size()*4);
    for(size_t i=0; i<pixels_in.size(); ++i)
    {
        for(size_t c=0;c<3;++c){
            pixels_out[4*i+c] = glm::uint8(std::pow(reinhard(pixels_in[i][c]),yamma)*255); 
        }
        pixels_out[4*i+3] = 255;
    }
    return pixels_out;
}