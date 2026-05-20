#include <stdint.h>
#include <bit>

class RNG {
private:
    uint32_t state;

public:
    RNG(uint32_t seed) : state(seed) {}

    uint32_t next_uint() {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        return state;
    }

    float next_float() {
        uint32_t u = (next_uint() >> 9) | 0x3f800000;
        return std::bit_cast<float>(u) - 1.0f;
    }
};