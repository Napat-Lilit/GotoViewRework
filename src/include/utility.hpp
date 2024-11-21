
#pragma once

#include <cmath>
#include <limits>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <memory>

constexpr float Pi = 3.14159265358979323846;
constexpr float InvPi = 0.31830988618379067154;
constexpr float Inv4Pi = 0.07957747154594766788;
constexpr float Sqrt2 = 1.41421356237309504880;
constexpr float Infinity = std::numeric_limits<float>::infinity();
constexpr float MachineEpsilon = std::numeric_limits<float>::epsilon() * 0.5;

// Better error handling, 
inline uint32_t FloatToBits(float f) {
    uint32_t ui;
    memcpy(&ui, &f, sizeof(float));
    return ui;
}
inline float BitsToFloat(uint32_t ui) {
    float f;
    memcpy(&f, &ui, sizeof(uint32_t));
    return f;
}

inline float NextFloatUp(float v) {
    // Handle infinity and negative zero for _NextFloatUp()_
    if (std::isinf(v) && v > 0.) return v;
    if (v == -0.f) v = 0.f;

    // Advance _v_ to next higher float
    uint32_t ui = FloatToBits(v);
    if (v >= 0)
        ++ui;
    else
        --ui;
    return BitsToFloat(ui);
}
inline float NextFloatDown(float v) {
    // Handle infinity and positive zero for _NextFloatDown()_
    if (std::isinf(v) && v < 0.) return v;
    if (v == 0.f) v = -0.f;
    uint32_t ui = FloatToBits(v);
    if (v > 0)
        --ui;
    else
        ++ui;
    return BitsToFloat(ui);
}

inline float Gamma(int n) {
    return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
}

inline float Clamp(float val, float low, float high) {
    if (val < low)
        return low;
    else if (val > high)
        return high;
    else
        return val;
}

inline float Radians(float deg) { return (Pi / 180.f) * deg; }
inline float Degrees(float rad) { return (180.f / Pi) * rad; }

// What you still lacks -> Literally random stuff
// 1. InverseRandomMax
// 2. random_float
// 3. random_int
