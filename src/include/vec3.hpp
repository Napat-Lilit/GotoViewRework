
#pragma once

#include <iostream>

#include "error.hpp"
#include "utility.hpp"

// Forward declaration of Normal3
class Normal3;

class Vec3 {
    public:
        Vec3() {x=y=z=0.f;}
        Vec3(float x, float y, float z) : x(x), y(y), z(z) {
            Assert(!HasNaNs());
        }
        bool HasNaNs() const {
            return std::isnan(x) || std::isnan(y) || std::isnan(z);
        }

        float operator[](int i) const { 
            Assert(i >= 0 && i <= 2);
            if (i == 0) return x;
            if (i == 1) return y;
            return z;
        }
        float &operator[](int i) { 
           Assert(i >= 0 && i <= 2);
           if (i == 0) return x;
           if (i == 1) return y;
           return z;
        }

        Vec3 operator+(const Vec3 &v) const {
            return Vec3(x + v.x, y + v.y, z + v.z);
        }
        Vec3& operator+=(const Vec3 &v) {
            x += v.x; y += v.y; z += v.z;
            return *this;
        }
        Vec3 operator-(const Vec3 &v) const {
            return Vec3(x - v.x, y - v.y, z - v.z);
        }
        Vec3& operator-=(const Vec3 &v) {
            x -= v.x; y -= v.y; z -= v.z;
            return *this;
        }
        bool operator==(const Vec3 &v) const {
            return x == v.x && y == v.y && z == v.z;
        }
        bool operator!=(const Vec3 &v) const {
            return x != v.x || y != v.y || z != v.z;
        }
        Vec3 operator*(float s) const { return Vec3(s*x, s*y, s*z); }
        Vec3 &operator*=(float s) {
            x *= s; y *= s; z *= s;
            return *this;
        }
        Vec3 operator/(float f) const {
            Assert(f != 0);
            float inv = (float)1 / f;
            return Vec3(x * inv, y * inv, z * inv);
        }
        Vec3 &operator/=(float f) {
            Assert(f != 0);
            float inv = (float)1 / f;
            x *= inv; y *= inv; z *= inv;
            return *this;
        }

        Vec3 operator-() const { return Vec3(-x, -y, -z); }
        float LengthSquared() const { return x * x + y * y + z * z; }
        float Length() const { return std::sqrt(LengthSquared()); }

        // explicit Vec3(const Normal3 &n) : x(n.x), y(n.y), z(n.z) {
        //     Assert(!n.HasNaNs());
        // }
        explicit Vec3(const Normal3 &n);    // Actual implementation is in cpp file

        float x,y,z;
};

// Type aliasing
// All of these are the exact same thing, namings are only for readability purpose
using Point3 = Vec3;
using RgbColor = Vec3;

// Normal must be handled seperately since transformation works a differently with normal
class Normal3 {
public:
    Normal3() { x = y = z = 0; }
    Normal3(float xx, float yy, float zz) : x(xx), y(yy), z(zz) {
        Assert(!HasNaNs());
    }
    bool HasNaNs() const {
        return std::isnan(x) || std::isnan(y) || std::isnan(z);
    }

    float operator[](int i) const { 
        Assert(i >= 0 && i <= 2);
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }
    
    float &operator[](int i) { 
        Assert(i >= 0 && i <= 2);
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }

    Normal3 operator+(const Normal3 &n) const {
        return Normal3(x + n.x, y + n.y, z + n.z);
    }
    Normal3& operator+=(const Normal3 &n) {
        x += n.x; y += n.y; z += n.z;
        return *this;
    }
    Normal3 operator- (const Normal3 &n) const {
        return Normal3(x - n.x, y - n.y, z - n.z);
    }
    Normal3& operator-=(const Normal3 &n) {
        x -= n.x; y -= n.y; z -= n.z;
        return *this;
    }
    bool operator==(const Normal3 &n) const {
        return x == n.x && y == n.y && z == n.z;
    }
    bool operator!=(const Normal3 &n) const {
        return x != n.x || y != n.y || z != n.z;
    }
    Normal3 operator*(float f) const {
        return Normal3(f*x, f*y, f*z);
    }
    Normal3 &operator*=(float f) {
        x *= f; y *= f; z *= f;
        return *this;
    }
    Normal3 operator/(float f) const {
        Assert(f != 0);
        float inv = (float)1 / f;
        return Normal3(x * inv, y * inv, z * inv);
    }
    Normal3 &operator/=(float f) {
        Assert(f != 0);
        float inv = (float)1 / f;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }

    Normal3 operator-() const {
        return Normal3(-x, -y, -z);
    }
    float LengthSquared() const { return x*x + y*y + z*z; }
    float Length() const        { return std::sqrt(LengthSquared()); }

    explicit Normal3(const Vec3 &v) : x(v.x), y(v.y), z(v.z) {
        Assert(!v.HasNaNs());
    }

    float x, y, z;
};

// Utility function
inline std::ostream& operator<<(std::ostream& out, const Vec3& v){
    out << "[" << v.x << ", " << v.y << ", " << v.z << "]";
    return out;
}
inline std::ostream& operator<<(std::ostream& out, const Normal3& n){
    out << "[" << n.x << ", " << n.y << ", " << n.z << "]";
    return out;
}

inline Vec3 operator*(float s, const Vec3 &v) {
    return v * s;
}
inline Normal3 operator*(float f, const Normal3 &n) {
    return Normal3(f * n.x, f * n.y, f * n.z);
}

// Maybe just inlining it will be clearner, but I digress
Vec3 Abs(const Vec3 &v);
Normal3 Abs(const Normal3 &v);

inline float Dot(const Normal3 &n1, const Vec3 &v2) {
    Assert(!n1.HasNaNs() && !v2.HasNaNs());
    return n1.x * v2.x + n1.y * v2.y + n1.z * v2.z;
}
inline float Dot(const Vec3 &v1, const Normal3 &n2) {
    Assert(!v1.HasNaNs() && !n2.HasNaNs());
    return v1.x * n2.x + v1.y * n2.y + v1.z * n2.z;
}
inline float Dot(const Normal3 &n1, const Normal3 &n2) {
    Assert(!n1.HasNaNs() && !n2.HasNaNs());
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}
inline float AbsDot(const Normal3 &n1, const Vec3 &v2) {
    Assert(!n1.HasNaNs() && !v2.HasNaNs());
    return std::abs(n1.x * v2.x + n1.y * v2.y + n1.z * v2.z);
}
inline float AbsDot(const Vec3 &v1, const Normal3 &n2) {
    Assert(!v1.HasNaNs() && !n2.HasNaNs());
    return std::abs(v1.x * n2.x + v1.y * n2.y + v1.z * n2.z);
}
inline float AbsDot(const Normal3 &n1, const Normal3 &n2) {
    Assert(!n1.HasNaNs() && !n2.HasNaNs());
    return std::abs(n1.x * n2.x + n1.y * n2.y + n1.z * n2.z);
}

inline Vec3 Cross(const Vec3 &v1, const Vec3 &v2) {
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vec3((v1y * v2z) - (v1z * v2y), 
                (v1z * v2x) - (v1x * v2z),
                (v1x * v2y) - (v1y * v2x));
}
inline Vec3 Cross(const Normal3 &v1, const Vec3 &v2) {
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vec3((v1y * v2z) - (v1z * v2y), 
                (v1z * v2x) - (v1x * v2z),
                (v1x * v2y) - (v1y * v2x));
}
inline Vec3 Cross(const Vec3 &v1, const Normal3 &v2) {
    Assert(!v1.HasNaNs() && !v2.HasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vec3((v1y * v2z) - (v1z * v2y), 
                (v1z * v2x) - (v1x * v2z),
                (v1x * v2y) - (v1y * v2x));
}

inline Vec3 Normalize(const Vec3 &v) {
    return v / v.Length();
}
inline Normal3 Normalize(const Normal3 &n) {
    return n / n.Length();
}

float MinComponent(const Vec3 &v);
float MaxComponent(const Vec3 &v);
int MaxDimension(const Vec3 &v);
Vec3 Min(const Vec3 &p1, const Vec3 &p2);
Vec3 Max(const Vec3 &p1, const Vec3 &p2);
Vec3 Permute(const Vec3 &p, int x, int y, int z);

// For generating a coordiante system from normal
inline void CoordinateSystem(const Vec3 &v1, Vec3 *v2, Vec3 *v3) {
    if (std::abs(v1.x) > std::abs(v1.y))
        *v2 = Vec3(-v1.z, 0, v1.x) / std::sqrt(v1.x * v1.x + v1.z * v1.z);
    else
        *v2 = Vec3(0, v1.z, -v1.y) / std::sqrt(v1.y * v1.y + v1.z * v1.z);
    *v3 = Cross(v1, *v2);
}

inline Point3 OffsetRayOrigin(const Point3 &p, const Vec3 &pError, const Normal3 &n, const Vec3 &w) {
    float d = Dot(Abs(n), pError);
    Vec3 offset = d * Vec3(n);

    if (Dot(w, n) < 0) offset = -offset;
    Point3 po = p + offset;
    // Round offset point _po_ away from _p_
    for (int i = 0; i < 3; ++i) {
        if (offset[i] > 0)
            po[i] = NextFloatUp(po[i]);
        else if (offset[i] < 0)
            po[i] = NextFloatDown(po[i]);
    }
    return po;
}

// What you still lacks
// 1. random_in_unit_sphere
// 2. random_unit_vector
// 3. random_in_unit_disk
// 4. reflect and refract
// 5. Offset ray origin
