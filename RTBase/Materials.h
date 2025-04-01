#pragma once

#include "Core.h"
#include "Imaging.h"
#include "Sampling.h"

#pragma warning( disable : 4244)

class BSDF;

class ShadingData
{
public:
	Vec3 x;
	Vec3 wo;
	Vec3 sNormal;
	Vec3 gNormal;
	float tu;
	float tv;
	Frame frame;
	BSDF* bsdf;
	float t;
	ShadingData() {}
	ShadingData(Vec3 _x, Vec3 n)
	{
		x = _x;
		gNormal = n;
		sNormal = n;
		bsdf = NULL;
	}

};

class ShadingHelper
{
public:
	static float fresnelDielectric(float cosTheta, float iorInt, float iorExt)
	{
		// Determine if we're entering or exiting the medium
		bool entering = cosTheta > 0.0f;

		// Swap indices of refraction if we're exiting the medium
		if (!entering) {
			std::swap(iorInt, iorExt);
			cosTheta = std::abs(cosTheta);
		}

		float eta = iorExt / iorInt;

		// Calculate the critical condition for total internal reflection
		float sinThetaT2 = eta * eta * (1.0f - cosTheta * cosTheta);

		// Check for total internal reflection
		if (sinThetaT2 > 1.0f) {
			return 1.0f; // Total reflection
		}

		// Calculate transmitted angle cosine
		float cosThetaT = std::sqrt(1.0f - sinThetaT2);

		// Fresnel equations for unpolarized light (Schlick's approximation)
		float r0 = std::pow((iorExt - iorInt) / (iorExt + iorInt), 2.0f);

		// Compute reflection probability
		float reflectance = r0 + (1.0f - r0) * std::pow(1.0f - cosTheta, 5.0f);

		return reflectance;
	}
	static Colour fresnelConductor(float cosTheta, Colour ior, Colour k)
	{
		// Ensure cosTheta is non-negative
		cosTheta = std::abs(cosTheta);

		// Compute common trigonometric terms
		float cos2 = cosTheta * cosTheta;
		float sin2 = 1.0f - cos2;

		// Compute squared values for parallel and perpendicular components
		Colour parallel_sq, perp_sq;

		// Parallel component squared
		parallel_sq.r = (ior.r * ior.r + k.r * k.r) * cos2 - 2 * ior.r * cosTheta + sin2;
		parallel_sq.g = (ior.g * ior.g + k.g * k.g) * cos2 - 2 * ior.g * cosTheta + sin2;
		parallel_sq.b = (ior.b * ior.b + k.b * k.b) * cos2 - 2 * ior.b * cosTheta + sin2;

		parallel_sq = parallel_sq / ((ior.r * ior.r + k.r * k.r) * cos2 + 2 * ior.r * cosTheta + sin2);

		// Perpendicular component squared 
		perp_sq.r = (ior.r * ior.r + k.r * k.r - 2 * ior.r * cosTheta + cos2);
		perp_sq.g = (ior.g * ior.g + k.g * k.g - 2 * ior.g * cosTheta + cos2);
		perp_sq.b = (ior.b * ior.b + k.b * k.b - 2 * ior.b * cosTheta + cos2);

		perp_sq = perp_sq / ((ior.r * ior.r + k.r * k.r) + 2 * ior.r * cosTheta + cos2);

		// Average the squared values
		return (parallel_sq + perp_sq) * 0.5f;
	}
	struct Vec3 {
		float x, y, z;

		// Dot product method
		float dot(const Vec3& other) const {
			return x * other.x + y * other.y + z * other.z;
		}

		// Absolute value of dot product with itself (equivalent to cos²θ)
		float absDot() const {
			return dot(*this);
		}
	};

	static float lambdaGGX(const Vec3& wi, float alpha)
	{
		// Calculate the projected area in the given direction
		float absDot = wi.absDot();

		// Compute the transformed view direction
		float tan2Theta = (1.0f - absDot) / absDot;

		// Lambda calculation using the GGX/Smith geometric shadowing function
		return 0.5f * (-1.0f + sqrt(1.0f + (alpha * alpha) * tan2Theta));
	}
	static float Gggx(Vec3 wi, Vec3 wo, float alpha)
	{
		// Compute the Smith Masking-Shadowing function
		// Using the height-correlated version of Smith geometric shadowing
		float lambdaWi = lambdaGGX(wi, alpha);
		float lambdaWo = lambdaGGX(wo, alpha);

		// G(wi, wo) = 1 / (1 + lambda(wi) + lambda(wo))
		return 1.0f / (1.0f + lambdaWi + lambdaWo);
	}

	static float Dggx(Vec3 h, float alpha)
	{
		// GGX/Trowbridge-Reitz microfacet distribution
		float alpha2 = alpha * alpha;
		float cosTheta = h.dot(Vec3{ 0, 1, 0 }); // Assuming y is the normal
		float cosTheta2 = cosTheta * cosTheta;

		// D(h) = α² / (π * cos⁴(θ) * (α² + tan²(θ))²)
		float tan2Theta = (1.0f - cosTheta2) / cosTheta2;

		return (alpha2) /
			(M_PI * cosTheta2 *
				(alpha2 + tan2Theta) *
				(alpha2 + tan2Theta));
	}
};

class BSDF
{
public:
	Colour emission;
	virtual Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf) = 0;
	virtual Colour evaluate(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual float PDF(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual bool isPureSpecular() = 0;
	virtual bool isTwoSided() = 0;
	bool isLight()
	{
		return emission.Lum() > 0 ? true : false;
	}
	void addLight(Colour _emission)
	{
		emission = _emission;
	}
	Colour emit(const ShadingData& shadingData, const Vec3& wi)
	{
		return emission;
	}
	virtual float mask(const ShadingData& shadingData) = 0;
};


class DiffuseBSDF : public BSDF
{
public:
	Texture* albedo;
	DiffuseBSDF() = default;
	DiffuseBSDF(Texture* _albedo)
	{
		albedo = _albedo;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Add correct sampling code here
		Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		pdf = wi.z / M_PI;
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		wi = shadingData.frame.toWorld(wi);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Add correct evaluation code here
		return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Add correct PDF code here
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		return SamplingDistributions::cosineHemispherePDF(wiLocal);
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class MirrorBSDF : public BSDF
{
public:
	Texture* albedo;
	MirrorBSDF() = default;
	MirrorBSDF(Texture* _albedo)
	{
		albedo = _albedo;
	}

	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);

		Vec3 wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);

		Vec3 wi = shadingData.frame.toWorld(wiLocal);

		float cosTheta = wiLocal.z; 
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / cosTheta;

		pdf = 1.0f;

		return wi;
	}

	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		return Colour(0.0f, 0.0f, 0.0f); 
	}

	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{

		return 0.0f;
	}

	bool isPureSpecular()
	{
		return true;
	}

	bool isTwoSided()
	{
		return true;
	}

	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};


class ConductorBSDF : public BSDF
{
public:
	Texture* albedo;
	Colour eta;
	Colour k;
	float alpha;
	ConductorBSDF() = default;
	ConductorBSDF(Texture* _albedo, Colour _eta, Colour _k, float roughness)
	{
		albedo = _albedo;
		eta = _eta;
		k = _k;
		alpha = 1.62142f * sqrtf(roughness);
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Replace this with Conductor sampling code
		Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		pdf = wi.z / M_PI;
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		wi = shadingData.frame.toWorld(wi);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Conductor evaluation code
		return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Conductor PDF
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		return SamplingDistributions::cosineHemispherePDF(wiLocal);
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class GlassBSDF : public BSDF
{
public:
	Texture* albedo;
	float intIOR;
	float extIOR;
	GlassBSDF() = default;
	GlassBSDF(Texture* _albedo, float _intIOR, float _extIOR)
	{
		albedo = _albedo;
		intIOR = _intIOR;
		extIOR = _extIOR;
	}

	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& refractedColour, float& pdf)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		float eta = extIOR / intIOR;
		float cosTheta = woLocal.z;

		// Calculate refraction direction
		float etaRatio = (cosTheta > 0) ? eta : 1.0f / eta;
		float sinThetaTSquared = etaRatio * etaRatio * (1.0f - cosTheta * cosTheta);

		if (sinThetaTSquared >= 1.0f) {
			// Total internal reflection - treat as perfect mirror
			Vec3 wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
			refractedColour = albedo->sample(shadingData.tu, shadingData.tv);
			pdf = 1.0f;
			return shadingData.frame.toWorld(wiLocal);
		}

		float cosThetaT = sqrtf(1.0f - sinThetaTSquared);
		if (cosTheta > 0) cosThetaT = -cosThetaT;

		Vec3 wiLocal = Vec3(-etaRatio * woLocal.x, -etaRatio * woLocal.y, cosThetaT);

		// Manual normalization
		float length = sqrtf(wiLocal.x * wiLocal.x + wiLocal.y * wiLocal.y + wiLocal.z * wiLocal.z);
		wiLocal.x /= length;
		wiLocal.y /= length;
		wiLocal.z /= length;

		// Fresnel term for energy conservation
		float F = fresnelDielectric(fabsf(cosTheta), extIOR, intIOR);

		// Energy compensation for transmitted light
		float factor = (etaRatio * etaRatio) / (extIOR * extIOR) * (1.0f - F);
		refractedColour = albedo->sample(shadingData.tu, shadingData.tv) * factor;
		pdf = 1.0f;
		return shadingData.frame.toWorld(wiLocal);
	}

	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Glass is a delta material - evaluation is zero except for exact refraction direction
		return Colour(0.0f, 0.0f, 0.0f);
	}

	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Delta distribution - zero everywhere except exact refraction direction
		return 0.0f;
	}

	bool isPureSpecular() { return true; }
	bool isTwoSided() { return true; }

	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}

private:
	float fresnelDielectric(float cosThetaI, float etaI, float etaT)
	{
		// Manual clamp
		if (cosThetaI > 1.0f) cosThetaI = 1.0f;
		if (cosThetaI < -1.0f) cosThetaI = -1.0f;

		bool entering = cosThetaI > 0.0f;
		if (!entering) {
			std::swap(etaI, etaT);
			cosThetaI = fabsf(cosThetaI);
		}

		float sinThetaI = sqrtf(fmaxf(0.0f, 1.0f - cosThetaI * cosThetaI));
		float sinThetaT = etaI / etaT * sinThetaI;
		if (sinThetaT >= 1.0f) return 1.0f;

		float cosThetaT = sqrtf(fmaxf(0.0f, 1.0f - sinThetaT * sinThetaT));

		float Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
			((etaT * cosThetaI) + (etaI * cosThetaT));
		float Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
			((etaI * cosThetaI) + (etaT * cosThetaT));

		return (Rparl * Rparl + Rperp * Rperp) / 2.0f;
	}
};

class DielectricBSDF : public BSDF
{
public:
	Texture* albedo;
	float intIOR;
	float extIOR;
	float alpha;
	DielectricBSDF() = default;
	DielectricBSDF(Texture* _albedo, float _intIOR, float _extIOR, float roughness)
	{
		albedo = _albedo;
		intIOR = _intIOR;
		extIOR = _extIOR;
		alpha = 1.62142f * sqrtf(roughness);
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Replace this with Dielectric sampling code
		Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		pdf = wi.z / M_PI;
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		wi = shadingData.frame.toWorld(wi);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Dielectric evaluation code
		return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Dielectric PDF
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		return SamplingDistributions::cosineHemispherePDF(wiLocal);
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return false;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class OrenNayarBSDF : public BSDF
{
public:
	Texture* albedo;
	float sigma;

	OrenNayarBSDF() = default;
	OrenNayarBSDF(Texture* _albedo, float _sigma)
		: albedo(_albedo), sigma(_sigma) {
	}

	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());

		wi = shadingData.frame.toWorld(wi);

		pdf = wi.z / M_PI;
		reflectedColour = evaluateOrenNayar(shadingData, wi);

		return wi;
	}

	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		return evaluateOrenNayar(shadingData, wi);
	}

	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		return SamplingDistributions::cosineHemispherePDF(wiLocal);
	}

	bool isPureSpecular() { return false; }
	bool isTwoSided() { return true; }

	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}

private:
	Colour evaluateOrenNayar(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 wo = shadingData.frame.toLocal(-shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);

		float sigma2 = sigma * sigma;
		float A = 1.0f - (sigma2 / (2.0f * (sigma2 + 0.33f)));
		float B = (0.45f * sigma2) / (sigma2 + 0.09f);

		float cosTheta_i = std::max(0.0f, wiLocal.z);
		float cosTheta_o = std::max(0.0f, wo.z);

		float sinTheta_i = std::sqrt(1.0f - cosTheta_i * cosTheta_i);
		float sinTheta_o = std::sqrt(1.0f - cosTheta_o * cosTheta_o);

		Vec3 v1(wiLocal.x, wiLocal.y, 0);
		Vec3 v2(wo.x, wo.y, 0);

		float l1 = std::sqrt(v1.x * v1.x + v1.y * v1.y);
		float l2 = std::sqrt(v2.x * v2.x + v2.y * v2.y);

		v1 = (l1 > 0) ? Vec3(v1.x / l1, v1.y / l1, 0) : Vec3(0, 0, 0);
		v2 = (l2 > 0) ? Vec3(v2.x / l2, v2.y / l2, 0) : Vec3(0, 0, 0);

		float maxAngle = std::max(0.0f, v1.x * v2.x + v1.y * v2.y);
		float cosDeltaPhi = maxAngle;

		Colour albedoColor = albedo->sample(shadingData.tu, shadingData.tv);

		float reflectance = (A + B * cosDeltaPhi * sinTheta_i * sinTheta_o /
			std::max(cosTheta_i, cosTheta_o)) / M_PI;

		return albedoColor * reflectance;
	}
};

class PlasticBSDF : public BSDF
{
public:
	Texture* albedo;
	float intIOR;
	float extIOR;
	float alpha;
	PlasticBSDF() = default;
	PlasticBSDF(Texture* _albedo, float _intIOR, float _extIOR, float roughness)
	{
		albedo = _albedo;
		intIOR = _intIOR;
		extIOR = _extIOR;
		alpha = 1.62142f * sqrtf(roughness);
	}
	float alphaToPhongExponent()
	{
		return (2.0f / SQ(std::max(alpha, 0.001f))) - 2.0f;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Replace this with Plastic sampling code
		Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		pdf = wi.z / M_PI;
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		wi = shadingData.frame.toWorld(wi);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Plastic evaluation code
		return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Plastic PDF
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		return SamplingDistributions::cosineHemispherePDF(wiLocal);
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class LayeredBSDF : public BSDF
{
public:
	BSDF* base;
	Colour sigmaa;
	float thickness;
	float intIOR;
	float extIOR;
	LayeredBSDF() = default;
	LayeredBSDF(BSDF* _base, Colour _sigmaa, float _thickness, float _intIOR, float _extIOR)
	{
		base = _base;
		sigmaa = _sigmaa;
		thickness = _thickness;
		intIOR = _intIOR;
		extIOR = _extIOR;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Add code to include layered sampling
		return base->sample(shadingData, sampler, reflectedColour, pdf);
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Add code for evaluation of layer
		return base->evaluate(shadingData, wi);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Add code to include PDF for sampling layered BSDF
		return base->PDF(shadingData, wi);
	}
	bool isPureSpecular()
	{
		return base->isPureSpecular();
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return base->mask(shadingData);
	}
};