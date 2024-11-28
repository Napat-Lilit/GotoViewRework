
#include "vec3.hpp"

// This is only to prevent circular referencing between Vec3 and Normal3 explicit conversion
Vec3::Vec3(const Normal3 &n) : x(n.x), y(n.y), z(n.z) {
    Assert(!n.HasNaNs());    
}

Vec3 Abs(const Vec3 &v) {
    return Vec3(std::abs(v.x), std::abs(v.y), std::abs(v.z));
}
Normal3 Abs(const Normal3 &v) {
    return Normal3(std::abs(v.x), std::abs(v.y), std::abs(v.z));
}

float MinComponent(const Vec3 &v) {
    return std::min(v.x, std::min(v.y, v.z));
}
float MaxComponent(const Vec3 &v) {
    return std::max(v.x, std::max(v.y, v.z));
}
int MaxDimension(const Vec3 &v) {
    return (v.x > v.y) ? ((v.x > v.z) ? 0 : 2) : ((v.y > v.z) ? 1 : 2);
}
Vec3 Min(const Vec3 &p1, const Vec3 &p2) {
    return Vec3(std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z));
}
Vec3 Max(const Vec3 &p1, const Vec3 &p2) {
    return Vec3(std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z));
}

Vec3 Permute(const Vec3 &p, int x, int y, int z) {
    return Vec3(p[x], p[y], p[z]);
}

// --------------------------------------------------------------------------------------------------

Vec3 Sqrt(const Vec3 &s) {
    Vec3 ret;
    ret.x = std::sqrt(ret.x);
    ret.y = std::sqrt(ret.y);
    ret.z = std::sqrt(ret.z);
    Assert(!ret.HasNaNs());
    return ret;
}