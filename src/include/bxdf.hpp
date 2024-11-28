
#pragma once

#include <optional>

#include "vec3.hpp"

// 2 Points to keep in mind
// 1. We will not consider volume rendering here, like, at all
// 2. We don't care about hard screen here, namely, when 2 transparent objects intersect each other 

// Interesting notes:
// 1. Disney uber BSDF

// Specifying which type of surface this is -> Since we didn't do multiple importance sampling, it really necessary for us here ?
enum BxDFType {
    Unset = 0,
    BSDF_REFLECTION   = 1 << 0,
    BSDF_TRANSMISSION = 1 << 1,
    BSDF_DIFFUSE      = 1 << 2,
    BSDF_GLOSSY       = 1 << 3,
    BSDF_SPECULAR     = 1 << 4,
    BSDF_ALL          = BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR |
                        BSDF_REFLECTION | BSDF_TRANSMISSION,
};
struct BSDFSample {
    BSDFSample(BxDFType flags = Unset, RgbColor f = RgbColor(), Vec3 wi = Vec3(), 
    float pdf = 0.f, float eta = 1)
    : f(f), wi(wi), pdf(pdf), types(types), eta(eta) {}

    // Note that unset flag imply unsuccessful sampling
    BxDFType types;
    RgbColor f;
    Vec3 wi;
    float pdf;
    float eta;
};

// ------------------------------------------------------------- Utility functions

// Reflection Declarations
float FrDielectric(float cosThetaI, float etaI, float etaT);
RgbColor FrConductor(float cosThetaI, const RgbColor &etaI,
                     const RgbColor &etaT, const RgbColor &k);

// BSDF Inline Functions
inline float CosTheta(const Vec3 &w) { return w.z; }
inline float Cos2Theta(const Vec3 &w) { return w.z * w.z; }
inline float AbsCosTheta(const Vec3 &w) { return std::abs(w.z); }
inline float Sin2Theta(const Vec3 &w) {
    return std::max((float)0, (float)1 - Cos2Theta(w));
}
inline float SinTheta(const Vec3 &w) { return std::sqrt(Sin2Theta(w)); }
inline float TanTheta(const Vec3 &w) { return SinTheta(w) / CosTheta(w); }
inline float Tan2Theta(const Vec3 &w) {
    return Sin2Theta(w) / Cos2Theta(w);
}
inline float CosPhi(const Vec3 &w) {
    float sinTheta = SinTheta(w);
    return (sinTheta == 0) ? 1 : Clamp(w.x / sinTheta, -1, 1);
}
inline float SinPhi(const Vec3 &w) {
    float sinTheta = SinTheta(w);
    return (sinTheta == 0) ? 0 : Clamp(w.y / sinTheta, -1, 1);
}
inline float Cos2Phi(const Vec3 &w) { return CosPhi(w) * CosPhi(w); }
inline float Sin2Phi(const Vec3 &w) { return SinPhi(w) * SinPhi(w); }
inline float CosDPhi(const Vec3 &wa, const Vec3 &wb) {
    float waxy = wa.x * wa.x + wa.y * wa.y;
    float wbxy = wb.x * wb.x + wb.y * wb.y;
    if (waxy == 0 || wbxy == 0)
        return 1;
    return Clamp((wa.x * wb.x + wa.y * wb.y) / std::sqrt(waxy * wbxy), -1, 1);
}

inline Vec3 Reflect(const Vec3 &wo, const Vec3 &n) {
    return -wo + 2 * Dot(wo, n) * n;
}
inline bool Refract(const Vec3 &wi, const Normal3 &n, float eta,
                    Vec3 *wt) {
    // Compute $\cos \theta_\roman{t}$ using Snell's law
    float cosThetaI = Dot(n, wi);
    float sin2ThetaI = std::max(float(0), float(1 - cosThetaI * cosThetaI));
    float sin2ThetaT = eta * eta * sin2ThetaI;

    // Handle total internal reflection for transmission
    if (sin2ThetaT >= 1) return false;
    float cosThetaT = std::sqrt(1 - sin2ThetaT);
    *wt = eta * -wi + (eta * cosThetaI - cosThetaT) * Vec3(n);
    return true;
}

inline bool SameHemisphere(const Vec3 &w, const Vec3 &wp) {
    return w.z * wp.z > 0;
}
inline bool SameHemisphere(const Vec3 &w, const Normal3 &wp) {
    return w.z * wp.z > 0;
}

// --------------------------------------------------------------------------------------------- Fresnel classes

// Recall that fresnel returns a single value
class Fresnel {
    public:
        virtual ~Fresnel() {};
        virtual RgbColor Evaluate(float cosI) const = 0;
};

class FresnelConductor : public Fresnel {
public:
    // <<FresnelConductor Public Methods>> 
    RgbColor Evaluate(float cosThetaI) const;
    FresnelConductor(const RgbColor &etaI, const RgbColor &etaT,
        const RgbColor &k) : etaI(etaI), etaT(etaT), k(k) { }

    // ***************************************************************************************** You need a method to adjust etaI

private:
    RgbColor etaI, etaT, k;
};
class FresnelDielectric : public Fresnel {
public:
    // <<FresnelDielectric Public Methods>> 
    RgbColor Evaluate(float cosThetaI) const;
    FresnelDielectric(float etaI, float etaT) : etaI(etaI), etaT(etaT) { }

    // ***************************************************************************************** You need a method to adjust etaI

private:
    float etaI, etaT;
};

// ---------------------------------------------------------------------------------------------
// This BxDF guy will be in place of material of the old system

class BxDF {
    public:
        virtual ~BxDF() {}
        BxDF(BxDFType type) : type(type) { }
        bool MatchesFlags(BxDFType t) const {
            return (type & t) == type;
        }

        // RgbColor here in place of spectrum 
        virtual RgbColor f(const Vec3 &wo, const Vec3 &wi) const = 0;
        virtual BSDFSample Sample_f(const Vec3 &wo, SampleSet2 sample) const = 0;
        virtual float Pdf(const Vec3 &wi, const Vec3 &wo) const = 0;

        const BxDFType type;
};

class SpecularReflection : public BxDF {
    public:
        // This may not be the brightest of ideas but I am really afrid of using raw pointers here
        SpecularReflection(const RgbColor &R, std::shared_ptr<Fresnel> fresnel) 
            : BxDF(BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)), R(R),
                fresnel(fresnel) { }
        RgbColor f(const Vec3 &wo, const Vec3 &wi) const { 
            return RgbColor(0.f); 
        }
        // RgbColor Sample_f(const Vec3 &wo, Vec3 *wi, const Point2f &sample, Float *pdf, BxDFType *sampledType) const;
        BSDFSample Sample_f(const Vec3 &wo, SampleSet2 sample) const;
        float Pdf(const Vec3 &wo, const Vec3 &wi) const {
            return 0;
        }

    private:
        // <<SpecularReflection Private Data>> 
        const RgbColor R;
        const std::shared_ptr<Fresnel> fresnel; // unique_ptr may helps us a bit on performance, but better play it simple for now
};

class SpecularTransmission : public BxDF {
    public:
        // <<SpecularTransmission Public Methods>> 
        SpecularTransmission(const RgbColor &T, float etaA, float etaB) 
        : BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR)), T(T), etaA(etaA), etaB(etaB), fresnel(etaA, etaB) {}

        RgbColor f(const Vec3 &wo, const Vec3 &wi) const { 
            return RgbColor(0.f); 
        }
        // RgbColor Sample_f(const Vec3 &wo, Vec3 *wi, const Point2f &sample, float *pdf, BxDFType *sampledType) const;
        BSDFSample Sample_f(const Vec3 &wo, SampleSet2 sample) const;
        float Pdf(const Vec3 &wo, const Vec3 &wi) const {
            return 0;
        }

    private:
        // <<SpecularTransmission Private Data>> 
        // This T is actually pretty interesting
        const RgbColor T;
        const float etaA, etaB;
        const FresnelDielectric fresnel;
};

class FresnelSpecular : public BxDF {
    public:
        // <<FresnelSpecular Public Methods>> 
        FresnelSpecular(const RgbColor &R, const RgbColor &T, float etaA, float etaB) 
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_SPECULAR)),
        R(R), T(T), etaA(etaA), etaB(etaB), fresnel(etaA, etaB) { }
        
        RgbColor f(const Vec3 &wo, const Vec3 &wi) const { 
            return RgbColor(0.f); 
        }
        BSDFSample Sample_f(const Vec3 &wo, SampleSet2 sample) const;
        float Pdf(const Vec3 &wo, const Vec3 &wi) const {
            return 0;
        }

    private:
        // <<FresnelSpecular Private Data>> 
        const RgbColor R, T;
        const float etaA, etaB;
        const FresnelDielectric fresnel;
};

class LambertianReflection : public BxDF {
    public:
        // <<LambertianReflection Public Methods>> 
        LambertianReflection(const RgbColor &R) : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(R) { }
        RgbColor f(const Vec3 &wo, const Vec3 &wi) const;
        
        // **************************************************************************** I want these 2, implement themselve if needed be
        // BSDFSample Sample_f(const Vec3 &wo, SampleSet2 sample) const;
        // float Pdf(const Vec3 &wo, const Vec3 &wi) const

    private:
        const RgbColor R;
};

// // ---------------------------------------------------------------------------------------------

// // BSDF -> Just a small wrapper over BxDF, to transfer between shading and rendering space
// // The idea is, once we know from hit record which object got hit, we will gain imformation regarding normal and material
// // Using these 2 together, we will initialize a new BSDF, and decide what the new direction + spectrum return will be
// class BSDF {
//     public:
//         BSDF(Normal3 ns, std::shared_ptr<BxDF> bxdf) : bxdf(bxdf) {
//             // Build a shading frame
//             CoordinateSystem(ns, &x, &y);
//         } 

//     private:
//         std::shared_ptr<BxDF> bxdf;
//         Vec3 x,y,z;
// };
