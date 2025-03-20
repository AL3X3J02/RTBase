#pragma once

#include "Core.h"
#include <random>
#include <algorithm>
class Sampler
{
public:
	virtual float next() = 0;
};

class MTRandom : public Sampler
{
public:
	std::mt19937 generator;
	std::uniform_real_distribution<float> dist;
	MTRandom(unsigned int seed = 1) : dist(0.0f, 1.0f)
	{
		generator.seed(seed);
	}
	float next()
	{
		return dist(generator);
	}
};

// Note all of these distributions assume z-up coordinate system
class SamplingDistributions
{
public:
	static Vec3 uniformSampleHemisphere(float r1, float r2)
	{

		float theta = std::acos(r1);  
		float phi = 2.0f * M_PI * r2;  

		float sin_theta = std::sin(theta);
		float cos_theta = std::cos(theta);
		float sin_phi = std::sin(phi);
		float cos_phi = std::cos(phi);

		float x = sin_theta * cos_phi;
		float y = sin_theta * sin_phi;
		float z = cos_theta;

		return Vec3(x, y, z);
	}

	static float uniformHemispherePDF(const Vec3 wi)
	{
		if (wi.z < 0)
			return 0.0f;

		return 1.0f / (2.0f * M_PI);
	}
	static Vec3 cosineSampleHemisphere(float r1, float r2)
	{
		float theta = std::acos(std::sqrt(r1));  
		float phi = 2.0f * M_PI * r2;  

		float sin_theta = std::sin(theta);
		float cos_theta = std::cos(theta);
		float sin_phi = std::sin(phi);
		float cos_phi = std::cos(phi);

		float x = sin_theta * cos_phi;
		float y = sin_theta * sin_phi;
		float z = cos_theta;

		return Vec3(x, y, z);
	}

	static float cosineHemispherePDF(const Vec3 wi)
	{
		if (wi.z < 0)
			return 0.0f;

		return wi.z / M_PI; 
	}
	static Vec3 uniformSampleSphere(float r1, float r2)
	{
		float theta = std::acos(1.0f - 2.0f * r1);  
		float phi = 2.0f * M_PI * r2;  

		float sin_theta = std::sin(theta);
		float cos_theta = std::cos(theta);
		float sin_phi = std::sin(phi);
		float cos_phi = std::cos(phi);

		float x = sin_theta * cos_phi;
		float y = sin_theta * sin_phi;
		float z = cos_theta;

		return Vec3(x, y, z);
	}

	static float uniformSpherePDF(const Vec3& wi)
	{
		return 1.0f / (4.0f * M_PI);
	}
};