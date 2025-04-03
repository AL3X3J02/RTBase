#pragma once

#include "Core.h"
#include "Imaging.h"
#include "Sampling.h"
#include <algorithm> // Add this header to support std::max and std::min functions

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

		// Constructor to convert from global Vec3
		Vec3() : x(0), y(0), z(0) {}
		Vec3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
		Vec3(const ::Vec3& v) : x(v.x), y(v.y), z(v.z) {}

		// Dot product method
		float dot(const Vec3& other) const {
			return x * other.x + y * other.y + z * other.z;
		}

		// Absolute value of dot product with itself (equivalent to cos²θ)
		float absDot() const {
			return dot(*this);
		}
	};

	// Helper method: Convert from global Vec3 to ShadingHelper::Vec3
	static Vec3 toInternalVec3(const ::Vec3& v) {
		return Vec3(v.x, v.y, v.z);
	}

	static float lambdaGGX(const Vec3& wi, float alpha)
	{
		// Calculate the projected area in the given direction
		float absDot = wi.absDot();

		// Compute the transformed view direction
		float tan2Theta = (1.0f - absDot) / absDot;

		// Lambda calculation using the GGX/Smith geometric shadowing function
		return 0.5f * (-1.0f + std::sqrt(1.0f + (alpha * alpha) * tan2Theta));
	}

	// GGX Distribution function for global Vec3
	static float Dggx(const ::Vec3& h, float alpha)
	{
		return Dggx(toInternalVec3(h), alpha);
	}

	// GGX Distribution function implementation for internal Vec3
	static float Dggx(const Vec3& h, float alpha)
	{
		// GGX/Trowbridge-Reitz microfacet distribution
		float alpha2 = alpha * alpha;
		float cosTheta = h.z; // Use h.z as the cosine of the angle with the normal
		float cosTheta2 = cosTheta * cosTheta;

		if (cosTheta <= 0) return 0.0f; // Prevent division by zero

		// D(h) = α² / (π * cos⁴(θ) * (α² + tan²(θ))²)
		float tan2Theta = (1.0f - cosTheta2) / (cosTheta2 + 1e-7f); // Avoid division by zero

		return (alpha2) /
			(M_PI * cosTheta2 *
				std::pow(alpha2 + tan2Theta, 2.0f));
	}

	// Calculate Smith shadowing function - version accepting global Vec3
	static float Gggx(const ::Vec3& wi, const ::Vec3& wo, float alpha)
	{
		return Gggx(toInternalVec3(wi), toInternalVec3(wo), alpha);
	}

	// Calculate Smith shadowing function - version accepting internal Vec3
	static float Gggx(const Vec3& wi, const Vec3& wo, float alpha)
	{
		return G1(wi, alpha) * G1(wo, alpha);
	}

private:
	// Lambda function for Smith G1 calculation
	static float Lambda(const Vec3& v, float alpha)
	{
		float cosTheta = std::abs(v.z);
		if (cosTheta <= 0) return 0.0f;
		
		float cos2Theta = cosTheta * cosTheta;
		float tan2Theta = (1.0f - cos2Theta) / cos2Theta;
		
		return 0.5f * (-1.0f + std::sqrt(1.0f + alpha * alpha * tan2Theta));
	}
	
	// Smith G1 function
	static float G1(const Vec3& v, float alpha)
	{
		return 1.0f / (1.0f + Lambda(v, alpha));
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

		// Perfect mirror reflection
		Vec3 wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
		wiLocal = wiLocal.normalize(); // Ensure direction vector is normalized

		Vec3 wi = shadingData.frame.toWorld(wiLocal);

		// Mirror BRDF contains the reciprocal of the cosine term, so no need to divide by cosTheta here
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv);

		pdf = 1.0f;
		return wi;
	}

	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Mirror is a Delta distribution, returns 0 in non-exact directions
		// But to avoid lighting calculation issues, return a very small value
		return albedo->sample(shadingData.tu, shadingData.tv) * 0.01f;
	}

	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// PDF for Delta distribution
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
	Colour eta;  // Complex refractive index of conductor (real part)
	Colour k;    // Extinction coefficient (imaginary part)
	float alpha; // Roughness
	ConductorBSDF() = default;
	ConductorBSDF(Texture* _albedo, Colour _eta, Colour _k, float roughness)
	{
		albedo = _albedo;
		eta = _eta;
		k = _k;
		alpha = 1.62142f * sqrtf(roughness); // Convert user input roughness to GGX alpha parameter
	}
	
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal;
		
		if (alpha < 0.01f) {
			// 光滑金属表面 - 使用镜面反射
			wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
			wiLocal = wiLocal.normalize();
		} else {
			// 粗糙金属表面 - 使用GGX微表面采样
			// 采样GGX分布生成微表面法线
			float r1 = sampler->next();
			float r2 = sampler->next();
			
			float phi = 2.0f * M_PI * r1;
			float cosTheta = sqrtf((1.0f - r2) / (1.0f + (alpha * alpha - 1.0f) * r2));
			float sinTheta = sqrtf(1.0f - cosTheta * cosTheta);
			
			// 微表面法线（局部空间）
			Vec3 m(sinTheta * cosf(phi), sinTheta * sinf(phi), cosTheta);
			
			// 根据微表面法线计算出射方向（通过反射）
			float dotWoM = woLocal.x * m.x + woLocal.y * m.y + woLocal.z * m.z;
			wiLocal = Vec3(2.0f * dotWoM * m.x - woLocal.x,
						   2.0f * dotWoM * m.y - woLocal.y,
						   2.0f * dotWoM * m.z - woLocal.z);
			
			// 确保在正确的半球
			if (wiLocal.z <= 0.0f) {
				// 反射到错误的半球，返回一个简单的镜面反射
				wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
				wiLocal = wiLocal.normalize();
			}
		}
		
		// 将局部方向转为世界方向
		Vec3 wi = shadingData.frame.toWorld(wiLocal);
		
		// 计算菲涅尔反射率，使用导体的菲涅尔公式
		float cosTheta = fabs(wiLocal.z);
		Colour F = ShadingHelper::fresnelConductor(cosTheta, eta, k);
		
		// 计算BRDF和PDF
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * F;
		pdf = 1.0f; // 简化：对于粗糙表面，这应该是GGX PDF，但这里简化处理
		
		return wi;
	}
	
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// 如果是几乎光滑的金属，我们认为是一个Delta分布
		if (alpha < 0.01f) {
			return albedo->sample(shadingData.tu, shadingData.tv) * 0.01f;
		}
		
		// 对于粗糙金属，使用微表面BRDF
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		
		// 计算半向量
		Vec3 h = (woLocal + wiLocal).normalize();
		
		// 计算法线分布函数D
		float D = ShadingHelper::Dggx(h, alpha);
		
		// 计算几何衰减G
		float G = ShadingHelper::Gggx(wiLocal, woLocal, alpha);
		
		// 计算菲涅尔项F
		float cosTheta = fabs(wiLocal.z);
		Colour F = ShadingHelper::fresnelConductor(cosTheta, eta, k);
		
		// 计算完整的微表面BRDF
		Colour f = albedo->sample(shadingData.tu, shadingData.tv) * D * G * F / (4.0f * woLocal.z * wiLocal.z);
		
		return f;
	}
	
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		if (alpha < 0.01f) {
			return 0.0f; // Delta分布
		}
		
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		
		// 计算半向量
		Vec3 h = (woLocal + wiLocal).normalize();
		
		// GGX分布PDF
		float D = ShadingHelper::Dggx(h, alpha);
		float Jacobian = 1.0f / (4.0f * std::abs(Dot(wiLocal, h)));
		
		return D * h.z * Jacobian;
	}
	
	bool isPureSpecular()
	{
		return alpha < 0.01f; // 仅当表面非常光滑时才认为是纯镜面
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
		float cosTheta = woLocal.z;
		bool entering = cosTheta > 0;
		
		// 确定介质之间的折射率比
		float etaI = entering ? extIOR : intIOR;
		float etaT = entering ? intIOR : extIOR;
		float eta = etaI / etaT;
		
		// 计算入射角的余弦和正弦值
		cosTheta = fabsf(cosTheta);
		float sinThetaI = sqrtf(fmaxf(0.0f, 1.0f - cosTheta * cosTheta));
		float sinThetaT = eta * sinThetaI;

		// 检查是否发生全内反射
		if (sinThetaT >= 1.0f) {
			// 全内反射 - 作为完美镜面处理
			Vec3 wiLocal = Vec3(-woLocal.x, -woLocal.y, entering ? woLocal.z : -woLocal.z);
			wiLocal = wiLocal.normalize();
			refractedColour = albedo->sample(shadingData.tu, shadingData.tv);
			pdf = 1.0f;
			return shadingData.frame.toWorld(wiLocal);
		}

		// 计算折射方向的余弦值
		float cosThetaT = sqrtf(fmaxf(0.0f, 1.0f - sinThetaT * sinThetaT));
		
		// 计算菲涅尔反射率
		float F = fresnelDielectric(cosTheta, etaI, etaT);
		
		// 随机决定是反射还是折射
		if (sampler->next() < F) {
			// 反射
			Vec3 wiLocal = Vec3(-woLocal.x, -woLocal.y, entering ? woLocal.z : -woLocal.z);
			wiLocal = wiLocal.normalize();
			refractedColour = albedo->sample(shadingData.tu, shadingData.tv);
			pdf = F;
			return shadingData.frame.toWorld(wiLocal);
		} else {
			// 折射
			// 计算折射方向
			Vec3 wiLocal;
			if (entering) {
				wiLocal = Vec3(-eta * woLocal.x, -eta * woLocal.y, -cosThetaT);
			} else {
				wiLocal = Vec3(-eta * woLocal.x, -eta * woLocal.y, cosThetaT);
			}
			wiLocal = wiLocal.normalize();
			
			// 计算能量补偿
			// 折射光的能量是(1-F)乘以折射率比的平方，用于确保能量守恒
			float energyScale = (1.0f - F);
			if (entering) {
				energyScale *= (etaT * etaT) / (etaI * etaI);
			}
			
			refractedColour = albedo->sample(shadingData.tu, shadingData.tv) * energyScale;
			pdf = 1.0f - F;
			return shadingData.frame.toWorld(wiLocal);
		}
	}

	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Glass is a delta material - evaluation is zero except for exact refraction direction
		// 由于直接光照计算可能会调用这个函数，我们返回一个非零但很小的值
		// 这不是物理上正确的，但可以防止玻璃显示为纯黑色
		return albedo->sample(shadingData.tu, shadingData.tv) * 0.01f;
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
	bool volumetricScattering;       // Whether to use volumetric scattering
	float scatteringCoefficient;     // Scattering coefficient
	
	DielectricBSDF() = default;
	DielectricBSDF(Texture* _albedo, float _intIOR, float _extIOR, float roughness)
	{
		albedo = _albedo;
		intIOR = _intIOR;
		extIOR = _extIOR;
		alpha = 1.62142f * sqrtf(roughness);
		volumetricScattering = false;  // 默认不启用体积散射
		scatteringCoefficient = 0.0f;  // 默认无散射
	}
	
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		float cosTheta = woLocal.z;
		bool entering = cosTheta > 0;
		
		// 确定介质之间的折射率比
		float etaI = entering ? extIOR : intIOR;
		float etaT = entering ? intIOR : extIOR;
		float eta = etaI / etaT;
		
		// 计算菲涅尔系数
		float cosThetaI = fabsf(cosTheta);
		float F = ShadingHelper::fresnelDielectric(cosThetaI, etaI, etaT);
		
		// 使用粗糙度采样微表面法线
		float r1 = sampler->next();
		float r2 = sampler->next();
		
		float phi = 2.0f * M_PI * r1;
		float roughness2 = alpha * alpha;
		float cosTheta2 = (1.0f - r2) / (1.0f + (roughness2 - 1.0f) * r2);
		float mCosTheta = sqrtf(cosTheta2);
		float mSinTheta = sqrtf(1.0f - cosTheta2);
		
		// 微表面法线（在局部空间中）
		Vec3 m(mSinTheta * cosf(phi), mSinTheta * sinf(phi), mCosTheta);
		
		// 决定是反射还是折射
		float dotWoM = woLocal.x * m.x + woLocal.y * m.y + woLocal.z * m.z;
		
		if (sampler->next() < F) {
			// 反射
			Vec3 wiLocal = Vec3(2.0f * dotWoM * m.x - woLocal.x,
							   2.0f * dotWoM * m.y - woLocal.y,
							   2.0f * dotWoM * m.z - woLocal.z);
			
			// 确保在正确的半球
			if (entering && wiLocal.z < 0) {
				// 如果反射到错误半球，使用镜面反射
				wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
				wiLocal = wiLocal.normalize();
			} else if (!entering && wiLocal.z > 0) {
				wiLocal = Vec3(-woLocal.x, -woLocal.y, -woLocal.z);
				wiLocal = wiLocal.normalize();
			}
			
			Vec3 wi = shadingData.frame.toWorld(wiLocal);
			
			// 计算反射的BRDF - 恢复简单的反射处理，防止反射变黑
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv);
			pdf = F;
			
			return wi;
		} else {
			// 折射
			// 计算折射方向
			float etaRatio = eta;
			float cosThetaT2 = 1.0f - etaRatio * etaRatio * (1.0f - dotWoM * dotWoM);
			
			if (cosThetaT2 <= 0.0f) {
				// 全内反射
				Vec3 wiLocal = Vec3(2.0f * dotWoM * m.x - woLocal.x,
								   2.0f * dotWoM * m.y - woLocal.y,
								   2.0f * dotWoM * m.z - woLocal.z);
				Vec3 wi = shadingData.frame.toWorld(wiLocal);
				reflectedColour = albedo->sample(shadingData.tu, shadingData.tv);
				pdf = 1.0f - F;
				return wi;
			}
			
			float cosThetaT = sqrtf(cosThetaT2);
			
			// 根据微表面法线计算折射方向
			Vec3 wiLocal;
			if (dotWoM > 0) {
				wiLocal = Vec3(-etaRatio * woLocal.x, -etaRatio * woLocal.y, -cosThetaT);
			} else {
				wiLocal = Vec3(-etaRatio * woLocal.x, -etaRatio * woLocal.y, cosThetaT);
			}
			wiLocal = wiLocal.normalize();
			
			Vec3 wi = shadingData.frame.toWorld(wiLocal);
			
			// --------------------------------------------------------------
			// 使用基于体积吸收原理的方法处理聚焦问题
			// --------------------------------------------------------------
			
			// 基础能量比例 - 保持较高透明度
			float energyScale = 5.0f; // 提高基础能量
			
			// 材质内部光线传播的物理效应
			if (entering) {
				energyScale *= (etaT * etaT) / (etaI * etaI);
			}
			
			// 检测是否处于聚焦区域
			float wiDotWo = fabsf(Dot(wiLocal.normalize(), woLocal.normalize()));
			
			// 模拟光线在介质中的体积吸收和散射
			float absorption = 0.15f; // 降低吸收率，提高透明度
			
			// 计算模拟光线在介质中传播的距离
			float travelDistance = 1.0f;
			if (wiDotWo > 0.7f) {
				// 在聚焦区域，光线传播距离更长，但增加得更缓和
				travelDistance += (wiDotWo - 0.7f) * 5.0f; // 从10.0降低到5.0
			}
			
			// 根据粗糙度调整体积散射
			float scattering = std::max(0.15f, alpha * 2.0f); // 提高最小散射值
			
			// 应用Beer-Lambert定律模拟体积吸收
			float transmittance = expf(-absorption * travelDistance / (scattering + 0.05f));
			
			// 在聚焦区域适度降低能量，但不要过度
			if (wiDotWo > 0.9f) { // 从0.85提高到0.9，缩小特殊处理区域
				float focusFactor = (wiDotWo - 0.9f) / 0.1f;
				transmittance *= std::max(0.3f, 1.0f - 0.5f * focusFactor); // 最多降低50%，最低保留30%
			}
			
			// 应用体积吸收效应
			energyScale *= transmittance;
			
			// 边缘处额外增强透明度，减少黑环
			if (cosThetaI < 0.3f) {
				energyScale *= 1.0f + (0.3f - cosThetaI) * 3.0f;
			}
			
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * energyScale;
			pdf = 1.0f - F;
			
			return wi;
		}
	}
	
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		
		bool entering = woLocal.z > 0;
		
		// 计算菲涅尔系数
		float etaI = entering ? extIOR : intIOR;
		float etaT = entering ? intIOR : extIOR;
		float cosThetaI = fabsf(woLocal.z);
		float F = ShadingHelper::fresnelDielectric(cosThetaI, etaI, etaT);
		
		// 对于粗糙的玻璃，使用微表面模型
		if (wiLocal.z * woLocal.z > 0) {
			// 同侧，可能是反射 - 返回简单的反射值，防止反射变黑
			if (alpha < 0.01f) {
				// 对于完全光滑的表面，返回一个小值避免黑点
				return albedo->sample(shadingData.tu, shadingData.tv) * 0.1f;
			}
			
			// 对于粗糙表面，计算正确的反射
			Vec3 h = (woLocal + wiLocal).normalize();
			float D = ShadingHelper::Dggx(h, alpha);
			float G = ShadingHelper::Gggx(wiLocal, woLocal, alpha);
			
			// 增加反射能量以避免黑色
			return albedo->sample(shadingData.tu, shadingData.tv) * F * D * G / (4.0f * fabsf(woLocal.z) * fabsf(wiLocal.z) + 1e-5f);
		} else {
			// 异侧，可能是折射
			// 确定介质之间的折射率比
			float eta = etaI / etaT;
			
			// 找到半向量（用于微表面模型）
			Vec3 htLocal;
			if (entering) {
				htLocal = -(wiLocal * eta + woLocal).normalize();
			} else {
				htLocal = -(woLocal * eta + wiLocal).normalize();
			}
			
			// 确保半向量指向正确方向
			if (htLocal.z < 0) htLocal = -htLocal;
			
			// 计算微表面分布和几何项
			float D = ShadingHelper::Dggx(htLocal, alpha);
			float G = ShadingHelper::Gggx(wiLocal, woLocal, alpha);
			
			// 计算折射的雅可比行列式
			float sqrtDenom = Dot(wiLocal, htLocal) + eta * Dot(woLocal, htLocal);
			float denom = std::max(sqrtDenom * sqrtDenom, 1e-5f);
			
			// 计算BTDF - 保持简化的折射计算
			float btdf = (1.0f - F) * D * G * eta * eta * fabsf(Dot(wiLocal, htLocal)) * fabsf(Dot(woLocal, htLocal));
			btdf /= std::max(fabsf(wiLocal.z * woLocal.z * denom), 1e-5f);
			
			// --------------------------------------------------------------
			// 使用基于体积吸收原理的方法处理聚焦问题
			// --------------------------------------------------------------
			
			// 基础能量比例 - 保持较高透明度
			float energyScale = 5.0f; // 提高基础能量
			
			// 材质内部光线传播的物理效应
			if (entering) {
				energyScale *= (etaT * etaT) / (etaI * etaI);
			}
			
			// 检测是否处于聚焦区域
			float wiDotWo = fabsf(Dot(wiLocal.normalize(), woLocal.normalize()));
			
			// 模拟光线在介质中的体积吸收和散射
			float absorption = 0.15f; // 降低吸收率，提高透明度
			
			// 计算模拟光线在介质中传播的距离
			float travelDistance = 1.0f;
			if (wiDotWo > 0.7f) {
				// 在聚焦区域，光线传播距离更长，但增加得更缓和
				travelDistance += (wiDotWo - 0.7f) * 5.0f; // 从10.0降低到5.0
			}
			
			// 根据粗糙度调整体积散射
			float scattering = std::max(0.15f, alpha * 2.0f); // 提高最小散射值
			
			// 应用Beer-Lambert定律模拟体积吸收
			float transmittance = expf(-absorption * travelDistance / (scattering + 0.05f));
			
			// 在聚焦区域适度降低能量，但不要过度
			if (wiDotWo > 0.9f) { // 从0.85提高到0.9，缩小特殊处理区域
				float focusFactor = (wiDotWo - 0.9f) / 0.1f;
				transmittance *= std::max(0.3f, 1.0f - 0.5f * focusFactor); // 最多降低50%，最低保留30%
			}
			
			// 应用体积吸收效应
			energyScale *= transmittance;
			
			// 边缘处额外增强透明度，减少黑环
			if (cosThetaI < 0.3f) {
				energyScale *= 1.0f + (0.3f - cosThetaI) * 3.0f;
			}
			
			return albedo->sample(shadingData.tu, shadingData.tv) * btdf * energyScale;
		}
	}
	
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		
		bool entering = woLocal.z > 0;
		
		// 计算菲涅尔系数
		float etaI = entering ? extIOR : intIOR;
		float etaT = entering ? intIOR : extIOR;
		float cosThetaI = fabsf(woLocal.z);
		float F = ShadingHelper::fresnelDielectric(cosThetaI, etaI, etaT);
		
		if (wiLocal.z * woLocal.z > 0) {
			// 反射
			return F;
		} else {
			// 折射
			return 1.0f - F;
		}
	}
	
	bool isPureSpecular()
	{
		return alpha < 0.01f;
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

class OrenNayarBSDF : public BSDF
{
public:
	Texture* albedo;
	float sigma;

	OrenNayarBSDF() = default;
	OrenNayarBSDF(Texture* _albedo, float _sigma)
		: albedo(_albedo), sigma(_sigma) {
		// 默认情况下，如果sigma传入null或0，使用一个合理的默认值
		if (sigma <= 0.0f) {
			sigma = 0.3f; // 轻微粗糙度的默认值
		}
	}

	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// 使用余弦加权采样选择出射方向
		Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		wi = shadingData.frame.toWorld(wi);

		// 余弦加权采样的PDF
		float cosTheta = Max(Dot(wi, shadingData.sNormal), 0.0f);
		pdf = cosTheta / M_PI;
		
		// 计算Oren-Nayar BRDF
		reflectedColour = evaluateOrenNayar(shadingData, wi);

		return wi;
	}

	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		return evaluateOrenNayar(shadingData, wi);
	}

	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// 余弦加权采样的PDF
		float cosTheta = Max(Dot(wi, shadingData.sNormal), 0.0f);
		return cosTheta / M_PI;
	}

	bool isPureSpecular() { return false; }
	bool isTwoSided() { return true; }

	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}

private:
	// 改进的Oren-Nayar模型实现
	Colour evaluateOrenNayar(const ShadingData& shadingData, const Vec3& wi)
	{
		// 将方向转换到局部坐标系
		Vec3 woLocal = shadingData.frame.toLocal(-shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		
		// 确保在上半球
		if (woLocal.z <= 0.0f || wiLocal.z <= 0.0f) {
			return Colour(0.0f, 0.0f, 0.0f);
		}

		// 计算Oren-Nayar模型参数
		float sigma2 = sigma * sigma;
		float A = 1.0f - (sigma2 / (2.0f * (sigma2 + 0.33f)));
		float B = (0.45f * sigma2) / (sigma2 + 0.09f);

		// 计算视角和光源角度的余弦值
		float cosTheta_i = wiLocal.z;  // 入射光方向与法线的夹角余弦
		float cosTheta_o = woLocal.z;  // 出射光方向与法线的夹角余弦
		
		// 计算正弦值
		float sinTheta_i = sqrtf(1.0f - cosTheta_i * cosTheta_i);
		float sinTheta_o = sqrtf(1.0f - cosTheta_o * cosTheta_o);
		
		// 计算入射光和出射光在切平面上的投影
		float cosPhiDiff = 0.0f;
		if (sinTheta_i > 1e-4f && sinTheta_o > 1e-4f) {
			// 归一化投影向量
			float invSinThetaI = 1.0f / sinTheta_i;
			float invSinThetaO = 1.0f / sinTheta_o;
			
			float cos_phi_i = wiLocal.x * invSinThetaI;
			float sin_phi_i = wiLocal.y * invSinThetaI;
			float cos_phi_o = woLocal.x * invSinThetaO;
			float sin_phi_o = woLocal.y * invSinThetaO;
			
			// 计算方位角差的余弦值
			cosPhiDiff = Max(0.0f, cos_phi_i * cos_phi_o + sin_phi_i * sin_phi_o);
		}
		
		// 计算alpha和beta
		float alpha = Max(cosTheta_i, cosTheta_o);
		float beta = Min(cosTheta_i, cosTheta_o);
		
		// 计算Oren-Nayar BRDF
		float L1 = A;
		float L2 = B * Max(0.0f, cosPhiDiff) * sinTheta_i * sinTheta_o / Max(cosTheta_i, cosTheta_o);
		
		// 最终BRDF值
		Colour albedoColor = albedo->sample(shadingData.tu, shadingData.tv);
		return albedoColor * (L1 + L2) / M_PI;
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
	
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		
		// 计算菲涅尔项，决定是反射还是漫反射
		float cosThetaI = fabsf(woLocal.z);
		float F = ShadingHelper::fresnelDielectric(cosThetaI, extIOR, intIOR);
		
		// 使用重要性采样选择漫反射或镜面反射
		if (sampler->next() < F) {
			// 镜面反射
			if (alpha < 0.01f) {
				// 完美镜面反射（光滑塑料）
				Vec3 wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
				wiLocal = wiLocal.normalize();
				Vec3 wi = shadingData.frame.toWorld(wiLocal);
				reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * F;
				pdf = F;
				return wi;
			} else {
				// 粗糙镜面反射（使用GGX微表面模型）
				// 采样微表面法线
				float r1 = sampler->next();
				float r2 = sampler->next();
				
				float phi = 2.0f * M_PI * r1;
				float cosTheta = sqrtf((1.0f - r2) / (1.0f + (alpha * alpha - 1.0f) * r2));
				float sinTheta = sqrtf(1.0f - cosTheta * cosTheta);
				
				// 微表面法线
				Vec3 m(sinTheta * cosf(phi), sinTheta * sinf(phi), cosTheta);
				
				// 计算反射方向
				float dotWoM = woLocal.x * m.x + woLocal.y * m.y + woLocal.z * m.z;
				Vec3 wiLocal = Vec3(2.0f * dotWoM * m.x - woLocal.x,
								   2.0f * dotWoM * m.y - woLocal.y,
								   2.0f * dotWoM * m.z - woLocal.z);
				
				// 确保在正确的半球
				if (wiLocal.z <= 0.0f) {
					// 反射到错误的半球，使用余弦采样
					wiLocal = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
				}
				
				Vec3 wi = shadingData.frame.toWorld(wiLocal);
				
				// 计算BRDF和PDF
				float D = ShadingHelper::Dggx(m, alpha);
				float G = ShadingHelper::Gggx(wiLocal, woLocal, alpha);
				
				reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * F * G * D / (4.0f * woLocal.z * wiLocal.z);
				
				// 计算正确的微表面PDF
				float jacobian = 1.0f / (4.0f * std::abs(Dot(wiLocal, m)));
				pdf = D * m.z * jacobian * F; // 考虑菲涅尔项影响的微表面PDF
				
				return wi;
			}
		} else {
			// 漫反射
			Vec3 wiLocal = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
			Vec3 wi = shadingData.frame.toWorld(wiLocal);
			
			// 缩放反照率以考虑能量守恒
			Colour diffuseAlbedo = albedo->sample(shadingData.tu, shadingData.tv) * (1.0f - F);
			reflectedColour = diffuseAlbedo / M_PI;
			
			pdf = wiLocal.z / M_PI * (1.0f - F);
			return wi;
		}
	}
	
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		
		// 计算菲涅尔系数
		float cosThetaI = fabsf(woLocal.z);
		float F = ShadingHelper::fresnelDielectric(cosThetaI, extIOR, intIOR);
		
		// 漫反射部分
		Colour diffuse = albedo->sample(shadingData.tu, shadingData.tv) * (1.0f - F) / M_PI;
		
		// 镜面反射部分
		Colour specular(0.0f, 0.0f, 0.0f);
		if (alpha < 0.01f) {
			// 光滑表面对漫反射光没有太大贡献，直接返回漫反射部分
			return diffuse;
		} else {
			// 计算镜面反射
			Vec3 h = (woLocal + wiLocal).normalize(); // 半向量
			
			float D = ShadingHelper::Dggx(h, alpha);
			float G = ShadingHelper::Gggx(wiLocal, woLocal, alpha);
			
			specular = albedo->sample(shadingData.tu, shadingData.tv) * F * D * G / (4.0f * woLocal.z * wiLocal.z);
		}
		
		return diffuse + specular;
	}
	
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		
		// 计算菲涅尔系数
		float cosThetaI = fabsf(woLocal.z);
		float F = ShadingHelper::fresnelDielectric(cosThetaI, extIOR, intIOR);
		
		if (alpha < 0.01f) {
			// 光滑表面，PDF是漫反射和镜面反射PDF的混合
			float diffusePDF = SamplingDistributions::cosineHemispherePDF(wiLocal) * (1.0f - F);
			// 镜面反射的PDF对于delta分布是0
			return diffusePDF;
		} else {
			// 粗糙表面，PDF是漫反射和微表面镜面反射PDF的混合
			float diffusePDF = SamplingDistributions::cosineHemispherePDF(wiLocal) * (1.0f - F);
			
			// 计算微表面镜面反射的PDF
			Vec3 h = (woLocal + wiLocal).normalize();
			float D = ShadingHelper::Dggx(h, alpha);
			float Jacobian = 1.0f / (4.0f * std::abs(Dot(wiLocal, h)));
			float specularPDF = D * h.z * Jacobian * F;
			
			return diffusePDF + specularPDF;
		}
	}
	
	bool isPureSpecular()
	{
		return false; // 塑料材质包含漫反射和镜面反射，不是纯镜面
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
		// 计算菲涅尔反射
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		float cosTheta = fabsf(woLocal.z);
		float F = ShadingHelper::fresnelDielectric(cosTheta, extIOR, intIOR);
		
		// 随机决定是表面反射还是穿透到基层
		if (sampler->next() < F) {
			// 表面反射
			Vec3 wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
			wiLocal = wiLocal.normalize();
			Vec3 wi = shadingData.frame.toWorld(wiLocal);
			reflectedColour = Colour(1.0f, 1.0f, 1.0f) * F; // 直接反射的光只受菲涅尔影响
			pdf = F;
			return wi;
		} else {
			// 穿透到基层
			// 从基层采样方向
			Colour baseReflectance;
			float basePdf;
			Vec3 wi = base->sample(shadingData, sampler, baseReflectance, basePdf);
			
			// 计算光线在涂层中传播的路径长度
			Vec3 wiLocal = shadingData.frame.toLocal(wi);
			float cosTheta_i = fabsf(wiLocal.z);
			float dist = thickness / cosTheta_i + thickness / cosTheta;
			
			// 应用Beer-Lambert法则计算衰减
			Colour transmission = Colour(
				expf(-sigmaa.r * dist),
				expf(-sigmaa.g * dist),
				expf(-sigmaa.b * dist)
			);
			
			// 考虑穿出涂层时的菲涅尔项
			float F_out = ShadingHelper::fresnelDielectric(cosTheta_i, intIOR, extIOR);
			
			// 最终反射率考虑进出涂层的透射率和基层反射率
			reflectedColour = baseReflectance * transmission * (1.0f - F) * (1.0f - F_out);
			pdf = (1.0f - F) * basePdf * (1.0f - F_out);
			
			return wi;
		}
	}
	
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		
		// 计算入射和出射角的余弦
		float cosTheta_o = fabsf(woLocal.z);
		float cosTheta_i = fabsf(wiLocal.z);
		
		// 菲涅尔项
		float F_in = ShadingHelper::fresnelDielectric(cosTheta_o, extIOR, intIOR);
		float F_out = ShadingHelper::fresnelDielectric(cosTheta_i, intIOR, extIOR);
		
		// 表面反射分量
		Colour specular(0.0f, 0.0f, 0.0f);
		if (abs(wiLocal.z - woLocal.z) < 0.001f && 
			abs(wiLocal.x + woLocal.x) < 0.001f && 
			abs(wiLocal.y + woLocal.y) < 0.001f) {
			// 镜面反射方向
			specular = Colour(1.0f, 1.0f, 1.0f) * F_in;
		}
		
		// 穿透到基层的分量
		Colour baseReflectance = base->evaluate(shadingData, wi);
		
		// 计算涂层中的传播距离
		float dist = thickness / cosTheta_i + thickness / cosTheta_o;
		
		// 应用Beer-Lambert法则
		Colour transmission = Colour(
			expf(-sigmaa.r * dist),
			expf(-sigmaa.g * dist),
			expf(-sigmaa.b * dist)
		);
		
		// 基层反射经过涂层的衰减
		Colour diffuse = baseReflectance * transmission * (1.0f - F_in) * (1.0f - F_out);
		
		return specular + diffuse;
	}
	
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		
		float cosTheta_o = fabsf(woLocal.z);
		float cosTheta_i = fabsf(wiLocal.z);
		
		// 菲涅尔项
		float F_in = ShadingHelper::fresnelDielectric(cosTheta_o, extIOR, intIOR);
		float F_out = ShadingHelper::fresnelDielectric(cosTheta_i, intIOR, extIOR);
		
		// 表面反射的PDF
		float specularPdf = 0.0f;
		if (abs(wiLocal.z - woLocal.z) < 0.001f && 
			abs(wiLocal.x + woLocal.x) < 0.001f && 
			abs(wiLocal.y + woLocal.y) < 0.001f) {
			// 镜面反射方向
			specularPdf = F_in;
		}
		
		// 基层的PDF
		float basePdf = base->PDF(shadingData, wi);
		
		// 考虑穿透概率
		float diffusePdf = (1.0f - F_in) * basePdf * (1.0f - F_out);
		
		return specularPdf + diffusePdf;
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