#pragma once

#include "Core.h"
#include "Imaging.h"
#include "Sampling.h"

#pragma warning( disable : 4244)

// å‰å‘è²æ˜Ž
class BSDF;
class Frame;

// ç¢ºä¿ShadingDataå®Œæ•´å®šç¾©åœ¨æ­¤è™•
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
		// Fresnel equations for dielectric materials (glass, water, etc.)
		float etaI = iorExt;
		float etaT = iorInt;
		
		// ç¢ºä¿cosThetaç‚ºæ­£å€¼
		bool entering = cosTheta > 0.0f;
		if (!entering) {
			std::swap(etaI, etaT);
			cosTheta = fabsf(cosTheta);
		}
		
		// è¨ˆç®—sin^2(thetaT)ï¼Œä½¿ç”¨Snellå®šå¾‹
		float sinThetaI = sqrtf(std::max(0.0f, 1 - cosTheta * cosTheta));
		float sinThetaT = etaI / etaT * sinThetaI;
		
		// å…¨å…§åå°„æª¢æŸ¥
		if (sinThetaT >= 1.0f) return 1.0f;
		
		float cosThetaT = sqrtf(std::max(0.0f, 1 - sinThetaT * sinThetaT));
		
		// è¨ˆç®—è²æ¶…çˆ¾æ–¹ç¨‹çš„å¹³è¡Œå’Œåž‚ç›´åˆ†é‡
		float Rs = ((etaT * cosTheta) - (etaI * cosThetaT)) / 
				  ((etaT * cosTheta) + (etaI * cosThetaT));
		float Rp = ((etaI * cosTheta) - (etaT * cosThetaT)) / 
				  ((etaI * cosTheta) + (etaT * cosThetaT));
		
		// è¿”å›žå¹³å‡å€¼çš„å¹³æ–¹ï¼Œè²æ¶…çˆ¾åå°„çŽ‡
		return 0.5f * (Rs * Rs + Rp * Rp);
	}
	
	static Colour fresnelConductor(float cosTheta, Colour eta, Colour k)
	{
		// Fresnel equations for conducting materials (metals)
		cosTheta = std::min(std::max(cosTheta, 0.0f), 1.0f);
		float cosTheta2 = cosTheta * cosTheta;
		float sinTheta2 = 1.0f - cosTheta2;
		
		Colour eta2 = Colour(eta.r * eta.r, eta.g * eta.g, eta.b * eta.b);
		Colour k2 = Colour(k.r * k.r, k.g * k.g, k.b * k.b);
		
		Colour t0 = Colour(eta2.r - k2.r - sinTheta2, eta2.g - k2.g - sinTheta2, eta2.b - k2.b - sinTheta2);
		
		Colour t0Sq = Colour(t0.r * t0.r, t0.g * t0.g, t0.b * t0.b);
		Colour k2eta2 = Colour(k2.r * eta2.r * 4.0f, k2.g * eta2.g * 4.0f, k2.b * eta2.b * 4.0f);
		Colour sum = Colour(t0Sq.r + k2eta2.r, t0Sq.g + k2eta2.g, t0Sq.b + k2eta2.b);
		
		Colour a2b2 = Colour(sqrtf(sum.r), sqrtf(sum.g), sqrtf(sum.b));
		
		Colour a = Colour(sqrtf((a2b2.r + t0.r) * 0.5f), sqrtf((a2b2.g + t0.g) * 0.5f), sqrtf((a2b2.b + t0.b) * 0.5f));
		
		Colour t1 = Colour(a2b2.r + cosTheta2, a2b2.g + cosTheta2, a2b2.b + cosTheta2);
		Colour t2 = Colour(a.r * cosTheta * 2.0f, a.g * cosTheta * 2.0f, a.b * cosTheta * 2.0f);
		
		Colour Rs = Colour((t1.r - t2.r) / (t1.r + t2.r), (t1.g - t2.g) / (t1.g + t2.g), (t1.b - t2.b) / (t1.b + t2.b));
		
		Colour t3 = Colour(cosTheta2 * a2b2.r + sinTheta2 * sinTheta2, cosTheta2 * a2b2.g + sinTheta2 * sinTheta2, cosTheta2 * a2b2.b + sinTheta2 * sinTheta2);
		Colour t4 = Colour(t2.r * sinTheta2, t2.g * sinTheta2, t2.b * sinTheta2);
		
		Colour Rp = Colour(Rs.r * (t3.r - t4.r) / (t3.r + t4.r), Rs.g * (t3.g - t4.g) / (t3.g + t4.g), Rs.b * (t3.b - t4.b) / (t3.b + t4.b));
		
		// è¿”å›žå¹³å‡å€¼
		return Colour((Rp.r + Rs.r) * 0.5f, (Rp.g + Rs.g) * 0.5f, (Rp.b + Rs.b) * 0.5f);
	}
	
	static float lambdaGGX(Vec3 wi, float alpha)
	{
		// Smithé®è”½å‡½æ•¸çš„è¼”åŠ©å‡½æ•¸
		float alpha2 = alpha * alpha;
		float cosTheta = fabsf(wi.z);
		float tanTheta2 = (1 - cosTheta * cosTheta) / (cosTheta * cosTheta);
		return 0.5f * (-1.0f + sqrtf(1.0f + alpha2 * tanTheta2));
	}
	
	static float Gggx(Vec3 wi, Vec3 wo, float alpha)
	{
		// Smithå¹¾ä½•è¡°æ¸›å› å­ï¼Œè€ƒæ…®äº†è¦–è§’å’Œå…‰æºæ–¹å‘
		return 1.0f / (1.0f + lambdaGGX(wi, alpha) + lambdaGGX(wo, alpha));
	}
	
	static float Dggx(Vec3 h, float alpha)
	{
		// GGXæ³•ç·šåˆ†å¸ƒå‡½æ•¸
		if (h.z <= 0) return 0.0f;
		
		float alpha2 = alpha * alpha;
		float cosTheta2 = h.z * h.z;
		float tanTheta2 = (1 - cosTheta2) / cosTheta2;
		
		float denom = M_PI * cosTheta2 * cosTheta2 * (alpha2 + tanTheta2);
		return alpha2 / denom;
	}
};

class BSDF
{
public:
	Colour emission;
	virtual Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& indirect, float& pdf) = 0;
	virtual Colour evaluate(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual float PDF(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual bool isPureSpecular() = 0;
	virtual bool isTwoSided() = 0;
	virtual bool isDiffuse() { return false; }
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
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& indirect, float& pdf)
	{
		// ä½¿ç”¨ä½™å¼¦åŠ æƒåŠçƒé‡‡æ ·
		float r1 = sampler->next();
		float r2 = sampler->next();
		
		// ä½¿ç”¨ä½™å¼¦åŠçƒé‡‡æ ·èŽ·å–æœ¬åœ°åæ ‡ä¸­çš„æ–¹å‘
		Vec3 localWi = SamplingDistributions::cosineSampleHemisphere(r1, r2);
		
		// è®¡ç®—PDF
		pdf = SamplingDistributions::cosineHemispherePDF(localWi);
		
		// è®¡ç®—åå°„é¢œè‰²
		indirect = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		
		// å°†æœ¬åœ°åæ ‡è½¬æ¢åˆ°ä¸–ç•Œåæ ‡
		Vec3 wi = shadingData.frame.toWorld(localWi);
		
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// å°‡ä¸–ç•Œåæ¨™è½‰æ›ç‚ºå±€éƒ¨åæ¨™
		Vec3 localWi = shadingData.frame.toLocal(wi);
		
		// è¨ˆç®—ä½™å¼¦é …ï¼ˆç¢ºä¿éžè² ï¼‰
		float cosTheta = std::max(0.0f, localWi.z);
		
		// æ­£ç¢ºçš„Lambertian BRDFéœ€è¦ä¹˜ä»¥ä½™å¼¦é …
		return albedo->sample(shadingData.tu, shadingData.tv) * cosTheta / M_PI;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// å°†ä¸–ç•Œåæ ‡è½¬æ¢å›žæœ¬åœ°åæ ‡
		Vec3 localWi = shadingData.frame.toLocal(wi);
		
		// å¦‚æžœåœ¨èƒŒé¢ï¼ŒPDFä¸º0
		if (localWi.z <= 0.0f) return 0.0f;
		
		// è®¡ç®—ä½™å¼¦åŠçƒPDF
		return SamplingDistributions::cosineHemispherePDF(localWi);
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
	bool isDiffuse() override
	{
		return true;
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
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& indirect, float& pdf)
	{
		// è®¡ç®—å®Œç¾Žé•œé¢åå°„æ–¹å‘
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		Vec3 localWi = Vec3(-wo.x, -wo.y, wo.z); // é•œé¢åå°„
		
		// é•œé¢BSDFçš„PDFæ˜¯deltaåˆ†å¸ƒ
		pdf = 1.0f;
		
		// åå°„é¢œè‰²
		indirect = albedo->sample(shadingData.tu, shadingData.tv);
		
		// è½¬æ¢å›žä¸–ç•Œåæ ‡
		return shadingData.frame.toWorld(localWi);
	}
	
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// é•œé¢æ˜¯deltaåˆ†å¸ƒï¼Œevaluateæ€»æ˜¯è¿”å›ž0
		// å®žé™…è´¡çŒ®é€šè¿‡sampleå‡½æ•°è®¡ç®—
		return Colour(0.0f, 0.0f, 0.0f);
	}
	
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// é•œé¢æ˜¯deltaåˆ†å¸ƒï¼ŒPDFæ— æ³•ç›´æŽ¥è®¡ç®—
        // å¯¹äºŽdeltaåˆ†å¸ƒï¼ŒPDFåº”è¯¥æ˜¯æ— ç©·å¤§
        // ä½†åœ¨å®žé™…å®žçŽ°ä¸­ï¼Œæˆ‘ä»¬è¿”å›ž0ï¼Œå› ä¸ºè¿™ä¸ªPDFåªä¼šåœ¨MISä¸­ä½¿ç”¨
        // è€ŒMISä¼šè·³è¿‡deltaåˆ†å¸ƒ
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
		// å°‡è¦–ç·šæ–¹å‘è½‰æ›ç‚ºå±€éƒ¨åæ¨™
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		
		// æ ¹æ“šGGXåˆ†å¸ƒæŽ¡æ¨£å¾®è¡¨é¢æ³•ç·š
		float r1 = sampler->next();
		float r2 = sampler->next();
		
		// å¾žGGXåˆ†å¸ƒæŽ¡æ¨£å¾®è¡¨é¢æ³•ç·š
		float phi = 2.0f * M_PI * r1;
		float tanTheta2 = alpha * alpha * r2 / (1.0f - r2);
		float cosTheta = 1.0f / sqrtf(1.0f + tanTheta2);
		float sinTheta = sqrtf(1.0f - cosTheta * cosTheta);
		
		// è¨ˆç®—å¾®è¡¨é¢æ³•ç·š
		Vec3 m = Vec3(sinTheta * cosf(phi), sinTheta * sinf(phi), cosTheta);
		
		// æ ¹æ“šå¾®è¡¨é¢æ³•ç·šè¨ˆç®—åå°„æ–¹å‘
		float coeff = 2.0f * Dot(wo, m);
		Vec3 localWi = Vec3(coeff * m.x - wo.x, coeff * m.y - wo.y, coeff * m.z - wo.z);
		
		// ç¢ºä¿åå°„æ–¹å‘åœ¨åŒä¸€åŠçƒ
		if (localWi.z * wo.z <= 0.0f) {
			// å¾®è¡¨é¢åå°„ä½¿åå°„æ–¹å‘ç©¿é€è¡¨é¢ï¼Œé€™ä¸æ‡‰è©²ç™¼ç”Ÿ
			pdf = 0.0f;
			reflectedColour = Colour(0.0f, 0.0f, 0.0f);
			return shadingData.frame.toWorld(localWi);
		}
		
		// è¨ˆç®—NDFæ¦‚çŽ‡å¯†åº¦
		float D = ShadingHelper::Dggx(m, alpha);
		pdf = D * cosTheta / (4.0f * Dot(wo, m));
		
		// ç¢ºä¿PDFæ­£ç¢º
		if (pdf <= 0.0f || !std::isfinite(pdf)) {
			pdf = 0.0f;
			reflectedColour = Colour(0.0f, 0.0f, 0.0f);
			return shadingData.frame.toWorld(localWi);
		}
		
		// è¨ˆç®—BRDFå€¼
		float G = ShadingHelper::Gggx(wo, localWi, alpha);
		Colour F = ShadingHelper::fresnelConductor(Dot(wo, m), eta, k);
		
		// å®Œæ•´çš„å¾®è¡¨é¢BRDF
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * F * G / (4.0f * wo.z * localWi.z);
		
		return shadingData.frame.toWorld(localWi);
	}
	
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// å°‡æ–¹å‘è½‰æ›ç‚ºå±€éƒ¨åæ¨™
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		Vec3 localWi = shadingData.frame.toLocal(wi);
		
		// ç¢ºä¿åå°„æ–¹å‘åœ¨åŒä¸€åŠçƒ
		if (localWi.z * wo.z <= 0.0f) {
			return Colour(0.0f, 0.0f, 0.0f);
		}
		
		// è¨ˆç®—åŠå‘é‡
		Vec3 h = Vec3(wo.x + localWi.x, wo.y + localWi.y, wo.z + localWi.z);
		float length = sqrtf(h.x*h.x + h.y*h.y + h.z*h.z);
		h = Vec3(h.x / length, h.y / length, h.z / length); // æ­¸ä¸€åŒ–
		
		// è¨ˆç®—å¾®è¡¨é¢BRDFçš„å„å€‹éƒ¨åˆ†
		float D = ShadingHelper::Dggx(h, alpha);
		float G = ShadingHelper::Gggx(wo, localWi, alpha);
		Colour F = ShadingHelper::fresnelConductor(Dot(wo, h), eta, k);
		
		return albedo->sample(shadingData.tu, shadingData.tv) * F * (D * G / (4.0f * wo.z * localWi.z));
	}
	
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// å°‡æ–¹å‘è½‰æ›ç‚ºå±€éƒ¨åæ¨™
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		Vec3 localWi = shadingData.frame.toLocal(wi);
		
		// ç¢ºä¿åå°„æ–¹å‘åœ¨åŒä¸€åŠçƒ
		if (localWi.z * wo.z <= 0.0f) {
			return 0.0f;
		}
		
		// è¨ˆç®—åŠå‘é‡
		Vec3 h = Vec3(wo.x + localWi.x, wo.y + localWi.y, wo.z + localWi.z);
		float length = sqrtf(h.x*h.x + h.y*h.y + h.z*h.z);
		h = Vec3(h.x / length, h.y / length, h.z / length); // æ­¸ä¸€åŒ–
		
		// è¨ˆç®—PDF
		float D = ShadingHelper::Dggx(h, alpha);
		float pdf = D * h.z / (4.0f * Dot(wo, h));
		
		return pdf;
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
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// å°‡è¦–ç·šæ–¹å‘è½‰æ›ç‚ºå±€éƒ¨åæ¨™
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		
		// è¨ˆç®—è²æ¶…çˆ¾é …
		float cosThetaI = wo.z;
		float F = ShadingHelper::fresnelDielectric(cosThetaI, intIOR, extIOR);
		
		// æ ¹æ“šè²æ¶…çˆ¾é …æ±ºå®šæ˜¯åå°„é‚„æ˜¯æŠ˜å°„
		float r = sampler->next();
		
		if (r < F) {
			// åå°„
			Vec3 localWi = Vec3(-wo.x, -wo.y, wo.z); // é¡é¢åå°„
			pdf = F;
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv);
			return shadingData.frame.toWorld(localWi);
		} else {
			// æŠ˜å°„
			float etaI = extIOR;
			float etaT = intIOR;
			
			// å¦‚æžœå…‰ç·šå¾žå…§éƒ¨å°„å‡ºï¼Œäº¤æ›IOR
			bool entering = cosThetaI > 0.0f;
			if (!entering) {
				std::swap(etaI, etaT);
				cosThetaI = -cosThetaI;
			}
			
			float eta = etaI / etaT;
			float sinThetaI = sqrtf(std::max(0.0f, 1.0f - cosThetaI * cosThetaI));
			float sinThetaT = eta * sinThetaI;
			
			// æª¢æŸ¥å…¨å…§åå°„
			if (sinThetaT >= 1.0f) {
				// ç™¼ç”Ÿå…¨å…§åå°„ï¼Œæ”¹ç‚ºåå°„
				Vec3 localWi = Vec3(-wo.x, -wo.y, wo.z); // é¡é¢åå°„
				pdf = 1.0f;
				reflectedColour = albedo->sample(shadingData.tu, shadingData.tv);
				return shadingData.frame.toWorld(localWi);
			}
			
			float cosThetaT = sqrtf(std::max(0.0f, 1.0f - sinThetaT * sinThetaT));
			
			// è¨ˆç®—æŠ˜å°„æ–¹å‘
			Vec3 localWi;
			if (entering) {
				localWi = Vec3(-eta * wo.x, -eta * wo.y, -cosThetaT);
			} else {
				localWi = Vec3(-eta * wo.x, -eta * wo.y, cosThetaT);
			}
			
			pdf = 1.0f - F;
			
			// æ ¹æ“šèƒ½é‡å®ˆæ†ï¼ŒæŠ˜å°„æ™‚éœ€è¦ä¹˜ä»¥eta^2
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * (eta * eta);
			
			return shadingData.frame.toWorld(localWi);
		}
	}
	
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// ç”±æ–¼çŽ»ç’ƒæ˜¯ç´”é¡é¢æè³ªï¼Œevaluateç¸½æ˜¯è¿”å›ž0
		// å¯¦éš›è²¢ç»é€šéŽsampleå‡½æ•¸è¨ˆç®—
		return Colour(0.0f, 0.0f, 0.0f);
	}
	
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// çŽ»ç’ƒæ˜¯deltaåˆ†å¸ƒï¼ŒPDFç„¡æ³•ç›´æŽ¥è¨ˆç®—
		// åœ¨å¯¦éš›å¯¦ç¾ä¸­ï¼Œæˆ‘å€‘è¿”å›ž0ï¼Œé€™è¡¨ç¤ºåœ¨MISä¸­ä¸è€ƒæ…®å®ƒ
		return 0.0f;
	}
	
	bool isPureSpecular()
	{
		return true;
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
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& indirect, float& pdf)
	{
		// è½¬æ¢å…¥å°„æ–¹å‘åˆ°å±€éƒ¨åæ ‡
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		
		// é‡‡æ ·å¾®è¡¨é¢æ³•çº¿
		float r1 = sampler->next();
		float r2 = sampler->next();
		
		// é‡‡æ ·GGXåˆ†å¸ƒèŽ·å–å¾®è¡¨é¢æ³•çº¿
		float phi = 2.0f * M_PI * r1;
		float cosTheta = sqrtf((1.0f - r2) / (1.0f + (alpha*alpha - 1.0f) * r2));
		float sinTheta = sqrtf(1.0f - cosTheta*cosTheta);
		
		// å¾®è¡¨é¢æ³•çº¿
		Vec3 m = Vec3(sinTheta * cosf(phi), sinTheta * sinf(phi), cosTheta);
		
		// è®¡ç®—è²æ¶…å°”é¡¹
		float F = ShadingHelper::fresnelDielectric(fabsf(Dot(wo, m)), intIOR, extIOR);
		
		// è®¡ç®—GGXåˆ†å¸ƒDå€¼
		float D = ShadingHelper::Dggx(m, alpha);
		
		// éšæœºå†³å®šæ˜¯åå°„è¿˜æ˜¯æŠ˜å°„
		if (sampler->next() < F)
		{
			// åå°„
			Vec3 localWi = Vec3(2.0f * Dot(wo, m) * m.x - wo.x,
                              2.0f * Dot(wo, m) * m.y - wo.y,
                              2.0f * Dot(wo, m) * m.z - wo.z);
			
			// ç¡®ä¿åå°„åœ¨æ­£ç¡®çš„åŠçƒ
			if (localWi.z * wo.z <= 0.0f) {
				pdf = 0.0f;
				indirect = Colour(0.0f, 0.0f, 0.0f);
				return Vec3(0.0f, 0.0f, 1.0f);
			}
			
			// è®¡ç®—PDF
			float jacobian = 1.0f / (4.0f * Dot(wo, m));
			pdf = D * cosTheta * jacobian * F; // ä¹˜ä»¥é€‰æ‹©åå°„çš„æ¦‚çŽ‡
			
			indirect = albedo->sample(shadingData.tu, shadingData.tv);
			return shadingData.frame.toWorld(localWi);
		}
		else
		{
			// æŠ˜å°„
			float eta = (wo.z > 0) ? extIOR/intIOR : intIOR/extIOR;
			float cosThetaT;
			Vec3 localWi;
			
			// å°è¯•è®¡ç®—æŠ˜å°„æ–¹å‘
			bool canRefract = refract(wo, m, eta, localWi);
			if (!canRefract)
			{
				// å…¨å†…åå°„ï¼Œå½“ä½œåå°„å¤„ç†
				localWi = Vec3(2.0f * Dot(wo, m) * m.x - wo.x,
                            2.0f * Dot(wo, m) * m.y - wo.y,
                            2.0f * Dot(wo, m) * m.z - wo.z);
		pdf = 1.0f;
			}
			else
			{
				// æŠ˜å°„æˆåŠŸ
				pdf = D * cosTheta * (1.0f - F); // PDFä¹˜ä»¥é€‰æ‹©æŠ˜å°„çš„æ¦‚çŽ‡
			}
			
			indirect = albedo->sample(shadingData.tu, shadingData.tv);
			return shadingData.frame.toWorld(localWi);
		}
	}
	
	// å¸®åŠ©å‡½æ•°ï¼šè®¡ç®—æŠ˜å°„æ–¹å‘
	bool refract(const Vec3& wi, const Vec3& n, float eta, Vec3& wt)
	{
		float cosThetaI = Dot(wi, n);
		float sin2ThetaI = std::max(0.0f, 1.0f - cosThetaI*cosThetaI);
		float sin2ThetaT = eta * eta * sin2ThetaI;
		
		// å…¨å†…åå°„æ£€æŸ¥
		if (sin2ThetaT >= 1.0f) return false;
		
		float cosThetaT = sqrtf(1.0f - sin2ThetaT);
		wt = Vec3(eta * (wi.x - n.x * cosThetaI) - n.x * cosThetaT,
               eta * (wi.y - n.y * cosThetaI) - n.y * cosThetaT,
               eta * (wi.z - n.z * cosThetaI) - n.z * cosThetaT);
		
		return true;
	}
	
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// ä¸ŽGGX + è²æ¶…å°”åå°„/æŠ˜å°„ç›¸å…³çš„BSDFè¯„ä¼°
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		Vec3 localWi = shadingData.frame.toLocal(wi);
		
		// æ£€æŸ¥æ˜¯åå°„è¿˜æ˜¯æŠ˜å°„
		bool isReflection = wo.z * localWi.z > 0.0f;
		
		if (isReflection)
		{
			// è®¡ç®—åŠç¨‹å‘é‡
			Vec3 h = Vec3(wo.x + localWi.x, wo.y + localWi.y, wo.z + localWi.z);
			float length = sqrtf(h.x*h.x + h.y*h.y + h.z*h.z);
			h = Vec3(h.x / length, h.y / length, h.z / length); // å½’ä¸€åŒ–
			
			// è®¡ç®—BSDF
			float D = ShadingHelper::Dggx(h, alpha);
			float G = ShadingHelper::Gggx(wo, localWi, alpha);
			float F = ShadingHelper::fresnelDielectric(fabsf(Dot(wo, h)), intIOR, extIOR);
			
			return albedo->sample(shadingData.tu, shadingData.tv) * F * (D * G / (4.0f * fabsf(wo.z * localWi.z)));
		}
		else
		{
			// æŠ˜å°„æƒ…å†µæ›´å¤æ‚ï¼Œè¿™é‡Œç®€åŒ–å¤„ç†
			return Colour(0.0f, 0.0f, 0.0f);
		}
	}
	
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// ä¸ŽGGX + è²æ¶…å°”åå°„/æŠ˜å°„ç›¸å…³çš„PDF
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		Vec3 localWi = shadingData.frame.toLocal(wi);
		
		// æ£€æŸ¥æ˜¯åå°„è¿˜æ˜¯æŠ˜å°„
		bool isReflection = wo.z * localWi.z > 0.0f;
		
		if (isReflection)
		{
			// è®¡ç®—åŠç¨‹å‘é‡
			Vec3 h = Vec3(wo.x + localWi.x, wo.y + localWi.y, wo.z + localWi.z);
			float length = sqrtf(h.x*h.x + h.y*h.y + h.z*h.z);
			h = Vec3(h.x / length, h.y / length, h.z / length); // å½’ä¸€åŒ–
			
			// è®¡ç®—åå°„PDF
			float cosTheta = h.z;
			float D = ShadingHelper::Dggx(h, alpha);
			float jacobian = 1.0f / (4.0f * Dot(wo, h));
			float F = ShadingHelper::fresnelDielectric(fabsf(Dot(wo, h)), intIOR, extIOR);
			
			return D * cosTheta * jacobian * F; // ä¹˜ä»¥é€‰æ‹©åå°„çš„æ¦‚çŽ‡
		}
		else
		{
			// æŠ˜å°„æƒ…å†µæ›´å¤æ‚ï¼Œè¿™é‡Œç®€åŒ–å¤„ç†
			return 0.0f;
		}
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
	{
		albedo = _albedo;
		sigma = _sigma;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& indirect, float& pdf)
	{
		// ä½¿ç”¨ä¸Žæ¼«åå°„ç›¸åŒçš„ä½™å¼¦åŠ æƒåŠçƒé‡‡æ ·
		float r1 = sampler->next();
		float r2 = sampler->next();
		
		// ä½¿ç”¨ä½™å¼¦åŠçƒé‡‡æ ·èŽ·å–æœ¬åœ°åæ ‡ä¸­çš„æ–¹å‘
		Vec3 localWi = SamplingDistributions::cosineSampleHemisphere(r1, r2);
		
		// è®¡ç®—PDF
		pdf = SamplingDistributions::cosineHemispherePDF(localWi);
		
		// è®¡ç®—åå°„é¢œè‰²ï¼ˆOren-Nayaræ¨¡åž‹ï¼‰
		indirect = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		
		// å°†æœ¬åœ°åæ ‡è½¬æ¢åˆ°ä¸–ç•Œåæ ‡
		Vec3 wi = shadingData.frame.toWorld(localWi);
		
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// å°†æ–¹å‘è½¬æ¢åˆ°å±€éƒ¨åæ ‡
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		Vec3 localWi = shadingData.frame.toLocal(wi);
		
		// ç¡®ä¿åœ¨æ­£ç¡®çš„åŠçƒ
		if (localWi.z <= 0.0f) return Colour(0.0f, 0.0f, 0.0f);
		
		// è®¡ç®—Oren-Nayaråå°„æ¨¡åž‹
		float sigma2 = sigma * sigma;
		float A = 1.0f - 0.5f * (sigma2 / (sigma2 + 0.33f));
		float B = 0.45f * (sigma2 / (sigma2 + 0.09f));
		
		// è®¡ç®—è§’åº¦
		float thetai = acosf(localWi.z);
		float thetao = acosf(wo.z);
		
		// è®¡ç®—æ–¹ä½è§’å·®
		float cosphi = 0.0f;
		if (thetai > 0.0001f && thetao > 0.0001f) {
			cosphi = (localWi.x * wo.x + localWi.y * wo.y) / 
					(sqrtf(localWi.x * localWi.x + localWi.y * localWi.y) * 
					 sqrtf(wo.x * wo.x + wo.y * wo.y));
			cosphi = std::max(-1.0f, std::min(1.0f, cosphi));
		}
		
		// è®¡ç®—alphaå’Œbeta
		float alpha = std::max(thetai, thetao);
		float beta = std::min(thetai, thetao);
		
		// Oren-Nayaråå°„æ–¹ç¨‹
		float oren_nayar = A + B * std::max(0.0f, cosphi) * sinf(alpha) * tanf(beta);
		
		// æœ€ç»ˆåå°„é¢œè‰²
		return albedo->sample(shadingData.tu, shadingData.tv) * (oren_nayar / M_PI);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// OrenNayarä½¿ç”¨ä¸Žæ¼«åå°„ç›¸åŒçš„é‡‡æ ·ç­–ç•¥ï¼Œæ‰€ä»¥PDFä¹Ÿç›¸åŒ
		Vec3 localWi = shadingData.frame.toLocal(wi);
		
		// å¦‚æžœåœ¨èƒŒé¢ï¼ŒPDFä¸º0
		if (localWi.z <= 0.0f) return 0.0f;
		
		// è®¡ç®—ä½™å¼¦åŠçƒPDF
		return SamplingDistributions::cosineHemispherePDF(localWi);
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
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& indirect, float& pdf)
	{
		// å¡‘æ–™æè´¨åŒ…å«æ¼«åå°„å’Œé«˜å…‰ä¸¤ä¸ªæˆåˆ†
		// éœ€è¦å†³å®šé‡‡æ ·å“ªä¸ªæˆåˆ†
		
		// å‡è®¾æ¼«åå°„æ¯”ä¾‹ä¸º0.5ï¼Œé«˜å…‰æ¯”ä¾‹ä¸º0.5
		float diffuseRatio = 0.5f;
		
		if (sampler->next() < diffuseRatio)
		{
			// é‡‡æ ·æ¼«åå°„æˆåˆ†
			float r1 = sampler->next();
			float r2 = sampler->next();
			
			// ä½¿ç”¨ä½™å¼¦åŠçƒé‡‡æ ·
			Vec3 localWi = SamplingDistributions::cosineSampleHemisphere(r1, r2);
			pdf = SamplingDistributions::cosineHemispherePDF(localWi) * diffuseRatio;
			
			// è®¡ç®—æ¼«åå°„éƒ¨åˆ†
			indirect = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
			
			return shadingData.frame.toWorld(localWi);
		}
		else
		{
			// é‡‡æ ·é«˜å…‰æˆåˆ†
			Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
			
			// é‡‡æ ·å¾®è¡¨é¢æ³•çº¿
			float r1 = sampler->next();
			float r2 = sampler->next();
			
			// é‡‡æ ·GGXåˆ†å¸ƒ
			float phi = 2.0f * M_PI * r1;
			float cosTheta = sqrtf((1.0f - r2) / (1.0f + (alpha*alpha - 1.0f) * r2));
			float sinTheta = sqrtf(1.0f - cosTheta*cosTheta);
			
			Vec3 m = Vec3(sinTheta * cosf(phi), sinTheta * sinf(phi), cosTheta);
			
			// è®¡ç®—åå°„æ–¹å‘
			Vec3 localWi = Vec3(2.0f * Dot(wo, m) * m.x - wo.x,
                              2.0f * Dot(wo, m) * m.y - wo.y,
                              2.0f * Dot(wo, m) * m.z - wo.z);
			
			// ç¡®ä¿åå°„åœ¨æ­£ç¡®çš„åŠçƒ
			if (localWi.z <= 0.0f) {
				pdf = 0.0f;
				indirect = Colour(0.0f, 0.0f, 0.0f);
				return Vec3(0.0f, 0.0f, 1.0f);
			}
			
			// è®¡ç®—é«˜å…‰PDF
			float D = ShadingHelper::Dggx(m, alpha);
			float jacobian = 1.0f / (4.0f * Dot(wo, m));
			pdf = D * cosTheta * jacobian * (1.0f - diffuseRatio);
			
			// è®¡ç®—è²æ¶…å°”é¡¹
			float F = ShadingHelper::fresnelDielectric(fabsf(Dot(wo, m)), intIOR, extIOR);
			
			// è®¡ç®—é«˜å…‰è´¡çŒ®
			float G = ShadingHelper::Gggx(wo, localWi, alpha);
			float result = F * D * G / (4.0f * wo.z * localWi.z);
			indirect = Colour(result, result, result); // ä½¿ç”¨ç›¸åŒçš„RGBå€¼
			
			return shadingData.frame.toWorld(localWi);
		}
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		Vec3 localWi = shadingData.frame.toLocal(wi);
		
		// ç¡®ä¿åœ¨ä¸ŠåŠçƒ
		if (localWi.z <= 0.0f) return Colour(0.0f, 0.0f, 0.0f);
		
		// è®¡ç®—æ¼«åå°„éƒ¨åˆ†
		Colour diffuse = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		
		// è®¡ç®—é«˜å…‰éƒ¨åˆ†
		// è®¡ç®—åŠç¨‹å‘é‡
		Vec3 h = Vec3(wo.x + localWi.x, wo.y + localWi.y, wo.z + localWi.z);
		float length = sqrtf(h.x*h.x + h.y*h.y + h.z*h.z);
		h = Vec3(h.x / length, h.y / length, h.z / length); // å½’ä¸€åŒ–
		
		float D = ShadingHelper::Dggx(h, alpha);
		float G = ShadingHelper::Gggx(wo, localWi, alpha);
		float F = ShadingHelper::fresnelDielectric(fabsf(Dot(wo, h)), intIOR, extIOR);
		
		float result = F * D * G / (4.0f * wo.z * localWi.z);
		Colour specular = Colour(result, result, result); // ä½¿ç”¨ç›¸åŒçš„RGBå€¼
		
		// è¿”å›žæ¼«åå°„å’Œé«˜å…‰çš„ç»„åˆ
		return diffuse + specular;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		Vec3 localWi = shadingData.frame.toLocal(wi);
		
		// ç¡®ä¿åœ¨ä¸ŠåŠçƒ
		if (localWi.z <= 0.0f) return 0.0f;
		
		// å‡è®¾æ¼«åå°„æ¯”ä¾‹ä¸º0.5ï¼Œé«˜å…‰æ¯”ä¾‹ä¸º0.5
		float diffuseRatio = 0.5f;
		
		// è®¡ç®—æ¼«åå°„PDF
		float diffusePdf = SamplingDistributions::cosineHemispherePDF(localWi) * diffuseRatio;
		
		// è®¡ç®—é«˜å…‰PDF
		// è®¡ç®—åŠç¨‹å‘é‡
		Vec3 h = Vec3(wo.x + localWi.x, wo.y + localWi.y, wo.z + localWi.z);
		float length = sqrtf(h.x*h.x + h.y*h.y + h.z*h.z);
		h = Vec3(h.x / length, h.y / length, h.z / length); // å½’ä¸€åŒ–
		
		float D = ShadingHelper::Dggx(h, alpha);
		float cosTheta = h.z;
		float jacobian = 1.0f / (4.0f * Dot(wo, h));
		float specularPdf = D * cosTheta * jacobian * (1.0f - diffuseRatio);
		
		// è¿”å›žä¸¤ç§PDFçš„ç»„åˆ
		return diffusePdf + specularPdf;
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
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& indirect, float& pdf)
	{
		return base->sample(shadingData, sampler, indirect, pdf);
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		return base->evaluate(shadingData, wi);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
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
