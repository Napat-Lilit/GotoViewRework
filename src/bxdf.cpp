
#include "bxdf.hpp"

float FrDielectric(float cosThetaI, float etaI, float etaT) {
    cosThetaI = Clamp(cosThetaI, -1, 1);
    // Potentially swap indices of refraction
    bool entering = cosThetaI > 0.f;
    if (!entering) {
        std::swap(etaI, etaT);
        cosThetaI = std::abs(cosThetaI);
    }

    // Compute _cosThetaT_ using Snell's law
    float sinThetaI = std::sqrt(std::max((float)0, 1 - cosThetaI * cosThetaI));
    float sinThetaT = etaI / etaT * sinThetaI;

    // Handle total internal reflection
    if (sinThetaT >= 1) return 1;
    float cosThetaT = std::sqrt(std::max((float)0, 1 - sinThetaT * sinThetaT));
    float Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
                  ((etaT * cosThetaI) + (etaI * cosThetaT));
    float Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
                  ((etaI * cosThetaI) + (etaT * cosThetaT));
    return (Rparl * Rparl + Rperp * Rperp) / 2;
}

// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
RgbColor FrConductor(float cosThetaI, const RgbColor &etai,
                     const RgbColor &etat, const RgbColor &k) {
    cosThetaI = Clamp(cosThetaI, -1, 1);
    RgbColor eta = etat / etai;
    RgbColor etak = k / etai;

    float cosThetaI2 = cosThetaI * cosThetaI;
    float sinThetaI2 = 1. - cosThetaI2;
    RgbColor eta2 = eta * eta;
    RgbColor etak2 = etak * etak;

    RgbColor t0 = eta2 - etak2 - RgbColor(sinThetaI2);
    RgbColor a2plusb2 = Sqrt(t0 * t0 + 4 * eta2 * etak2);
    RgbColor t1 = a2plusb2 + cosThetaI2;
    RgbColor a = Sqrt(0.5f * (a2plusb2 + t0));
    RgbColor t2 = (float)2 * cosThetaI * a;
    RgbColor Rs = (t1 - t2) / (t1 + t2);

    RgbColor t3 = cosThetaI2 * a2plusb2 + sinThetaI2 * sinThetaI2;
    RgbColor t4 = t2 * sinThetaI2;
    RgbColor Rp = Rs * (t3 - t4) / (t3 + t4);

    return 0.5 * (Rp + Rs);
}

// --------------------------------------------------------------------------------------------------------

RgbColor FresnelConductor::Evaluate(float cosThetaI) const {
    return FrConductor(std::abs(cosThetaI), etaI, etaT, k);
}
RgbColor FresnelDielectric::Evaluate(float cosThetaI) const {
    return FrDielectric(cosThetaI, etaI, etaT);
}

// --------------------------------------------------------------------------------------------------------

BSDFSample SpecularReflection::Sample_f(const Vec3 &wo, SampleSet2 sample) const {
    Vec3 wi = Vec3(-wo.x, -wo.y, wo.z);
    RgbColor f = fresnel->Evaluate(CosTheta(wi)) * R / AbsCosTheta(wi);
    float pdf = 1.f;
    // eta is not needed -> reflected back to what ever the ray comes from

    BSDFSample materialInteraction(type, f, wi, pdf);
    return materialInteraction;
}

BSDFSample SpecularTransmission::Sample_f(const Vec3 &wo, SampleSet2 sample) const {
    // <<Figure out which is incident and which is transmitted>> 
    bool entering = CosTheta(wo) > 0;
    float etaI = entering ? etaA : etaB;
    float etaT = entering ? etaB : etaA;

    // <<Compute ray direction for specular transmission>> 
    Vec3 wi;
    // Flipping normal to face the incoming wo
    Normal3 n(0, 0, 1);
    Normal3 FaceforwardNormal = Dot(n, wo) < 0.f ? -n : n;
    if (!Refract(wo, FaceforwardNormal, etaI / etaT, &wi)) {
        BSDFSample materialInteraction;

        // The default sample type is unset, thus unsuccessful sampling
        return materialInteraction;
    }

    float pdf = 1;
    RgbColor ft = T * (RgbColor(1.) - fresnel.Evaluate(CosTheta(wi)));
    // <<Account for non-symmetry with transmission to different medium>> 
    // **************************************************************************** Not 100 % sure here though, either to comment out or not
    // if (mode == TransportMode::Radiance)
    //    ft *= (etaI * etaI) / (etaT * etaT);

    RgbColor f = ft / AbsCosTheta(wi);

    BSDFSample materialInteraction(type, f, wi, pdf);
    return materialInteraction;
}

BSDFSample FresnelSpecular::Sample_f(const Vec3 &wo, SampleSet2 sample) const {

    float F = FrDielectric(CosTheta(wo), etaA, etaB);
    if (sample.Sample1 < F) {
        // Compute specular reflection for _FresnelSpecular_

        // Compute perfect specular reflection direction
        Vec3 wi(-wo.x, -wo.y, wo.z);
        BxDFType sampledType = BxDFType(BSDF_SPECULAR | BSDF_REFLECTION);
        float pdf = F;
        RgbColor f = F * R / AbsCosTheta(wi);
        return BSDFSample(sampledType, f, wi, pdf);

    } else {
        // Compute specular transmission for _FresnelSpecular_
        // Figure out which $\eta$ is incident and which is transmitted
        bool entering = CosTheta(wo) > 0;
        float etaI = entering ? etaA : etaB;
        float etaT = entering ? etaB : etaA;

        // Compute ray direction for specular transmission
        Vec3 wi;
        Normal3 n(0, 0, 1);
        Normal3 FaceforwardNormal = Dot(n, wo) < 0.f ? -n : n;
        if (!Refract(wo, FaceforwardNormal, etaI / etaT, &wi)) {
            BSDFSample materialInteraction;

            // The default sample type is unset, thus unsuccessful sampling
            return materialInteraction;
        }

        RgbColor ft = T * (1 - F);
        // Account for non-symmetry with transmission to different medium
        // **************************************************************************** Not 100 % sure here though, either to comment out or not
        // if (mode == TransportMode::Radiance)
        //     ft *= (etaI * etaI) / (etaT * etaT);
        BxDFType sampledType = BxDFType(BSDF_SPECULAR | BSDF_TRANSMISSION);
        float pdf = 1 - F;
        RgbColor f = ft / AbsCosTheta(wi);
        return BSDFSample(sampledType, f, wi, pdf);
    }
}

RgbColor LambertianReflection::f(const Vec3 &wo, const Vec3 &wi) const {
    return R * InvPi;
}
