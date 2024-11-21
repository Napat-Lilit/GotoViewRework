
#pragma once

#include "ray.hpp"

class Aabb {
    public:
        Aabb() {
            float minNum = std::numeric_limits<float>::lowest();
            float maxNum = std::numeric_limits<float>::max();
            pMin = Point3(maxNum, maxNum, maxNum);
            pMax = Point3(minNum, minNum, minNum);
        }
        explicit Aabb(const Point3 &p) : pMin(p), pMax(p) {}
        Aabb(const Point3 &p1, const Point3 &p2)
            : pMin(std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z)),
            pMax(std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z)) {}

        const Point3 &operator[](int i) const;
        Point3 &operator[](int i);

        bool operator==(const Aabb &b) const {
            return b.pMin == pMin && b.pMax == pMax;
        }
        bool operator!=(const Aabb &b) const {
            return b.pMin != pMin || b.pMax != pMax;
        }
        Vec3 Diagonal() const { return pMax - pMin; }
        float SurfaceArea() const {
            Vec3 d = Diagonal();
            return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
        }
        float Volume() const {
            Vec3 d = Diagonal();
            return d.x * d.y * d.z;
        }
        int MaximumExtent() const {
            Vec3 d = Diagonal();
            if (d.x > d.y && d.x > d.z)
                return 0;
            else if (d.y > d.z)
                return 1;
            else
                return 2;
        }

        bool IntersectP(const Ray &ray, float *hitt0 = nullptr, float *hitt1 = nullptr) const;
        // inline bool IntersectP(const Ray &ray, const Vec3 &invDir, const int dirIsNeg[3]) const;

        Point3 pMin, pMax;
};

inline const Point3 &Aabb::operator[](int i) const {
    Assert(i == 0 || i == 1);
    return (i == 0) ? pMin : pMax;
}
inline Point3 &Aabb::operator[](int i) {
    Assert(i == 0 || i == 1);
    return (i == 0) ? pMin : pMax;
}
Aabb Union(const Aabb &b, const Point3 &p) {
    Aabb ret;
    ret.pMin = Min(b.pMin, p);
    ret.pMax = Max(b.pMax, p);
    return ret;
}
Aabb Union(const Aabb &b1, const Aabb &b2) {
    Aabb ret;
    ret.pMin = Min(b1.pMin, b2.pMin);
    ret.pMax = Max(b1.pMax, b2.pMax);
    return ret;
}

// Note that this is IntersectP, not Intersect
inline bool Aabb::IntersectP(const Ray &ray, float *hitt0, float *hitt1) const {
    // float t0 = 0, t1 = ray.tMax;
    float t0 = 0, t1 = Infinity;

    for (int i = 0; i < 3; ++i) {
        // Update interval for _i_th bounding box slab
        float invRayDir = 1 / ray.dir[i];
        float tNear = (pMin[i] - ray.orig[i]) * invRayDir;
        float tFar = (pMax[i] - ray.orig[i]) * invRayDir;

        // Update parametric interval from slab intersection $t$ values
        if (tNear > tFar) std::swap(tNear, tFar);

        // Update _tFar_ to ensure robust ray--bounds intersection
        tFar *= 1 + 2 * Gamma(3);
        t0 = tNear > t0 ? tNear : t0;
        t1 = tFar < t1 ? tFar : t1;
        if (t0 > t1) return false;
    }
    if (hitt0) *hitt0 = t0;
    if (hitt1) *hitt1 = t1;
    return true;
}
