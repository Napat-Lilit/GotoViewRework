
#pragma once

#include "vec3.hpp"

class Ray {
  public:
    Ray(const Point3 &o, const Vec3 &d, const float IoR = 1.f) : orig(o), dir(d), currentIoR(IoR) {}

    Point3 operator()(float t) const { return orig + dir * t; }
    bool HasNaNs() const { return (orig.HasNaNs() || dir.HasNaNs()); }

    Point3 orig;
    Vec3 dir;
    // const Medium *medium;
    const float currentIoR;
};
