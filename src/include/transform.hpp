
#pragma once

#include "vec3.hpp"

struct Matrix4x4 {
    Matrix4x4() {
        m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.f;
        m[0][1] = m[0][2] = m[0][3] = m[1][0] =
        m[1][2] = m[1][3] = m[2][0] = m[2][1] = m[2][3] =
        m[3][0] = m[3][1] = m[3][2] = 0.f;
    }
    Matrix4x4(float mat[4][4]);
    Matrix4x4(float t00, float t01, float t02, float t03,
            float t10, float t11, float t12, float t13,
            float t20, float t21, float t22, float t23,
            float t30, float t31, float t32, float t33);

    static Matrix4x4 Mul(const Matrix4x4 &m1, const Matrix4x4 &m2) {
        Matrix4x4 r;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                r.m[i][j] = m1.m[i][0] * m2.m[0][j] + m1.m[i][1] * m2.m[1][j] +
                            m1.m[i][2] * m2.m[2][j] + m1.m[i][3] * m2.m[3][j];
        return r;
    }
    friend Matrix4x4 Transpose(const Matrix4x4 &m);
    friend Matrix4x4 Inverse(const Matrix4x4 &m);

    float m[4][4];
};

class Transform {
    public:
        Transform() { }
        Transform(const float mat[4][4]) {
            m = Matrix4x4(mat[0][0], mat[0][1], mat[0][2], mat[0][3],
                            mat[1][0], mat[1][1], mat[1][2], mat[1][3],
                            mat[2][0], mat[2][1], mat[2][2], mat[2][3],
                            mat[3][0], mat[3][1], mat[3][2], mat[3][3]);
            mInv = Inverse(m);
        }
        Transform(const Matrix4x4 &m) : m(m), mInv(Inverse(m)) { }
        Transform(const Matrix4x4 &m, const Matrix4x4 &mInv) 
            : m(m), mInv(mInv) {
        }

        Transform operator* (const Transform &t2) const;

        inline Point3 PointTransform(const Point3& p) const;
        inline Vec3 VecTransform(const Vec3& p) const;
        inline Normal3 NormalTransform(const Normal3& p) const;

    private:
        Matrix4x4 m, mInv;
};

// For composite transformation
Transform Transform :: operator*(const Transform &t2) const {
    return Transform(Matrix4x4::Mul(m, t2.m),
                     Matrix4x4::Mul(t2.mInv, mInv));
}

inline Point3 Transform::PointTransform(const Point3& p) const {
    float x = p.x, y = p.y, z = p.z;
    float xp = m.m[0][0]*x + m.m[0][1]*y + m.m[0][2]*z + m.m[0][3];
    float yp = m.m[1][0]*x + m.m[1][1]*y + m.m[1][2]*z + m.m[1][3];
    float zp = m.m[2][0]*x + m.m[2][1]*y + m.m[2][2]*z + m.m[2][3];
    float wp = m.m[3][0]*x + m.m[3][1]*y + m.m[3][2]*z + m.m[3][3];
    if (wp == 1) return Point3(xp, yp, zp);
    else         return Point3(xp, yp, zp) / wp;
}
inline Vec3 Transform::VecTransform(const Vec3& v) const {
    float x = v.x, y = v.y, z = v.z;
    return Vec3(m.m[0][0]*x + m.m[0][1]*y + m.m[0][2]*z,
                m.m[1][0]*x + m.m[1][1]*y + m.m[1][2]*z,
                m.m[2][0]*x + m.m[2][1]*y + m.m[2][2]*z);
}
inline Normal3 Transform::NormalTransform(const Normal3& n) const {
    float x = n.x, y = n.y, z = n.z;
    return Normal3(mInv.m[0][0]*x + mInv.m[1][0]*y + mInv.m[2][0]*z,
                mInv.m[0][1]*x + mInv.m[1][1]*y + mInv.m[2][1]*z,
                mInv.m[0][2]*x + mInv.m[1][2]*y + mInv.m[2][2]*z);
}

// Model transformations will happen before everything else
