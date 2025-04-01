#pragma once

#include "Core.h"
#include "Sampling.h"
#include "Geometry.h"
#include "Imaging.h"
#include "Materials.h"
#include "Lights.h"
#include "Scene.h"
#include "GamesEngineeringBase.h"
#include <thread>
#include <functional>
#include <OpenImageDenoise/oidn.hpp>  // Added OIDN header

static MTRandom* samplers;

class RayTracer
{
public:
	Scene* scene;
	GamesEngineeringBase::Window* canvas;
	Film* film;
	std::thread** threads;
	int numProcs;

	// Add buffers for denoiser
	float* colorBuffer;
	float* albedoBuffer;
	float* normalBuffer;
	float* outputBuffer;
	bool enableDenoising;

	// Denoiser strength adjuster
	// Move this into Main 
	float denoiseStrength = 0.0f; // Adjustable between 0.0 (no denoising) and 1.0 (full denoising). Unfortunatly this doesnt work too well making it way too blurry.

	void init(Scene* _scene, GamesEngineeringBase::Window* _canvas, bool _enableDenoising = true)
	{
		scene = _scene;
		canvas = _canvas;
		film = new Film();
		film->init((unsigned int)scene->camera.width, (unsigned int)scene->camera.height, new MitchellNetravaliFilter()); //Replace new with either BoxFilter, GaussianFilter or MitchellNetravaliFilter
		SYSTEM_INFO sysInfo;
		GetSystemInfo(&sysInfo);
		numProcs = sysInfo.dwNumberOfProcessors;
		threads = new std::thread * [numProcs];
		enableDenoising = _enableDenoising;

		// Allocate samplers array once
		samplers = new MTRandom[numProcs];
		// Initialize each sampler
		for (int i = 0; i < numProcs; i++) {
			samplers[i] = MTRandom();
		}

		// Allocate denoiser buffers
		int size = film->width * film->height * 3;
		colorBuffer = new float[size];
		albedoBuffer = new float[size];
		normalBuffer = new float[size];
		outputBuffer = new float[size];

		clear();
	}

	void clear()
	{
		film->clear();

		// Clear denoiser buffers
		int size = film->width * film->height * 3;
		memset(colorBuffer, 0, size * sizeof(float));
		memset(albedoBuffer, 0, size * sizeof(float));
		memset(normalBuffer, 0, size * sizeof(float));
		memset(outputBuffer, 0, size * sizeof(float));
	}

	Colour computeDirect(ShadingData shadingData, Sampler* sampler)
	{
		if (shadingData.bsdf->isPureSpecular() == true)
		{
			return Colour(0.0f, 0.0f, 0.0f);
		}
		// Sample a light
		float pmf;
		Light* light = scene->sampleLight(sampler, pmf);
		// Sample a point on the light
		float pdf;
		Colour emitted;
		Vec3 p = light->sample(shadingData, sampler, emitted, pdf);
		if (light->isArea())
		{
			// Calculate GTerm
			Vec3 wi = p - shadingData.x;
			float l = wi.lengthSq();
			wi = wi.normalize();
			float GTerm = (max(Dot(wi, shadingData.sNormal), 0.0f) * max(-Dot(wi, light->normal(shadingData, wi)), 0.0f)) / l;
			if (GTerm > 0)
			{
				// Trace
				if (scene->visible(shadingData.x, p))
				{
					// Shade
					return shadingData.bsdf->evaluate(shadingData, wi) * emitted * GTerm / (pmf * pdf);
				}
			}
		}
		else
		{
			// Calculate GTerm
			Vec3 wi = p;
			float GTerm = max(Dot(wi, shadingData.sNormal), 0.0f);
			if (GTerm > 0)
			{
				// Trace
				if (scene->visible(shadingData.x, shadingData.x + (p * 10000.0f)))
				{
					// Shade
					return shadingData.bsdf->evaluate(shadingData, wi) * emitted * GTerm / (pmf * pdf);
				}
			}
		}
		return Colour(0.0f, 0.0f, 0.0f);
	}

	Colour pathTrace(Ray& r, Colour& pathThroughput, int depth, Sampler* sampler, bool canHitLight = true)
	{
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		if (shadingData.t < FLT_MAX)
		{
			if (shadingData.bsdf->isLight())
			{
				if (canHitLight == true)
				{
					return pathThroughput * shadingData.bsdf->emit(shadingData, shadingData.wo);
				}
				else
				{
					return Colour(0.0f, 0.0f, 0.0f);
				}
			}
			Colour direct = pathThroughput * computeDirect(shadingData, sampler);
			if (depth > 5)
			{
				return direct;
			}
			float russianRouletteProbability = min(pathThroughput.Lum(), 0.9f);
			if (sampler->next() < russianRouletteProbability)
			{
				pathThroughput = pathThroughput / russianRouletteProbability;
			}
			else
			{
				return direct;
			}
			Colour indirect;
			float pdf;
			Vec3 wi = shadingData.bsdf->sample(shadingData, sampler, indirect, pdf);

			pathThroughput = pathThroughput * indirect * fabsf(Dot(wi, shadingData.sNormal)) / pdf;
			r.init(shadingData.x + (wi * EPSILON), wi);

			return (direct + pathTrace(r, pathThroughput, depth + 1, sampler, shadingData.bsdf->isPureSpecular()));
		}
		return scene->background->evaluate(shadingData, r.dir);
	}

	Colour direct(Ray& r, Sampler* sampler)
	{
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		if (shadingData.t < FLT_MAX)
		{
			if (shadingData.bsdf->isLight())
			{
				return shadingData.bsdf->emit(shadingData, shadingData.wo);
			}
			return computeDirect(shadingData, sampler);
		}
		return Colour(0.0f, 0.0f, 0.0f);
	}

	Colour albedo(Ray& r)
	{
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		if (shadingData.t < FLT_MAX)
		{
			if (shadingData.bsdf->isLight())
			{
				return shadingData.bsdf->emit(shadingData, shadingData.wo);
			}
			return shadingData.bsdf->evaluate(shadingData, Vec3(0, 1, 0));
		}
		return scene->background->evaluate(shadingData, r.dir);
	}

	Colour viewNormals(Ray& r)
	{
		IntersectionData intersection = scene->traverse(r);
		if (intersection.t < FLT_MAX)
		{
			ShadingData shadingData = scene->calculateShadingData(intersection, r);
			return Colour(fabsf(shadingData.sNormal.x), fabsf(shadingData.sNormal.y), fabsf(shadingData.sNormal.z));
		}
		return Colour(0.0f, 0.0f, 0.0f);
	}

	void connectToCamera(Vec3 p, Vec3 n, Colour col)
	{
		// Check if point p is visible from camera
		bool visible = scene->visible(p, scene->camera.origin);
		if (!visible)
			return;

		// Compute direction to camera and normalize it
		Vec3 wi = scene->camera.origin - p;
		float distanceSquared = wi.lengthSq();
		wi = wi.normalize();

		// Check if normal is facing the camera (cosine term)
		float cosTheta = max(Dot(n, wi), 0.0f);
		if (cosTheta <= 0.0f)
			return;

		// Project point onto camera - implement this or find equivalent in your camera class
		float x, y;
		// Create a simplified version of projectOntoCamera
		bool onCamera = projectOntoCamera(p, x, y);
		if (!onCamera)
			return;

		// Compute geometry term between p and camera
		float cameraCosTheta = max(Dot(-wi, scene->camera.viewDirection), 0.0f);

		// Geometry term: cos(θ) * cos(θ') / r²
		float G = (cosTheta * cameraCosTheta) / distanceSquared;

		// Compute weight for film splat
		float filmArea = scene->camera.width * scene->camera.height;
		float W_e = 1.0f / (filmArea * pow(cameraCosTheta, 4.0f));

		// Final contribution to splat onto film
		Colour contribution = col * G * W_e;

		// Splat onto film at projected coordinates
		film->splat(x, y, contribution);

		// Update display
		unsigned char r = (unsigned char)(contribution.r * 255);
		unsigned char g = (unsigned char)(contribution.g * 255);
		unsigned char b = (unsigned char)(contribution.b * 255);
		film->tonemap((unsigned int)x, (unsigned int)y, r, g, b);
		canvas->draw((unsigned int)x, (unsigned int)y, r, g, b);
	}

	// Helper method to project a world point onto camera plane
	bool projectOntoCamera(Vec3 p, float& x, float& y)
	{
		// Create a ray from point p to camera origin
		Vec3 dir = (scene->camera.origin - p).normalize();

		// Check if ray direction is facing the camera
		float dot = Dot(dir, scene->camera.viewDirection);
		if (dot >= 0) // Point is behind camera
			return false;

		// Calculate the distance to the camera plane
		float t = -Dot(p - scene->camera.origin, scene->camera.viewDirection) / dot;
		if (t <= 0) // Point is behind camera
			return false;

		// Calculate the intersection with camera plane
		Vec3 intersection = p + dir * t;

		// Convert to camera space
		Vec3 localPos = intersection - scene->camera.origin;

		// Project onto camera plane (assuming camera.up and camera.right are available)
		// If not available, you'll need to compute them based on viewDirection
		Vec3 right = Cross(scene->camera.viewDirection, Vec3(0, 1, 0)).normalize();
		Vec3 up = Cross(right, scene->camera.viewDirection).normalize();

		float u = Dot(localPos, right);
		float v = Dot(localPos, up);

		// Convert to normalized device coordinates
		float aspectRatio = scene->camera.width / scene->camera.height;
		float fov = 60.0f * (3.14159f / 180.0f); // Assuming 60 degree FOV, adjust as needed
		float scale = tan(fov * 0.5f);

		u = (u / scale / aspectRatio + 1.0f) * 0.5f * scene->camera.width;
		v = (-v / scale + 1.0f) * 0.5f * scene->camera.height;

		// Check if within screen bounds
		if (u < 0 || u >= scene->camera.width || v < 0 || v >= scene->camera.height)
			return false;

		x = u;
		y = v;
		return true;
	}

	void lightTrace(Sampler* sampler)
	{
		// Sample a light
		float pmf;
		Light* light = scene->sampleLight(sampler, pmf);
		if (!light)
			return;

		// Create a dummy ShadingData for light sampling
		ShadingData dummyShading;
		dummyShading.x = Vec3(0, 0, 0); // Position will be overwritten

		// Sample a point on the light (using the method from your existing code)
		float pdf;
		Colour emitted;
		Vec3 p = light->sample(dummyShading, sampler, emitted, pdf);

		// Get the light normal at this point
		Vec3 lightNormal;
		if (light->isArea()) {
			// For area lights, get the normal
			// Here we're creating a dummy direction toward the sampled point
			Vec3 dummyWi = Vec3(0, 1, 0); // Will be ignored for getting normal
			lightNormal = light->normal(dummyShading, dummyWi);
		}
		else {
			// For non-area lights, use a default up vector
			lightNormal = Vec3(0, 1, 0);
		}

		// Sample a direction from the light (hemisphere sampling around normal)
		Vec3 wi = sampleHemisphere(lightNormal, sampler);

		// Compute emission in this direction
		Colour Le = emitted / (pmf * pdf);

		// Create a ray starting at p in direction wi
		Ray r;
		r.init(p + (wi * EPSILON), wi);

		// Connect the initial point to the camera
		connectToCamera(p, lightNormal, Le);

		// Continue the light path
		Colour pathThroughput = Le;
		lightTracePath(r, pathThroughput, Le, sampler);
	}

	// Helper function for hemisphere sampling
	Vec3 sampleHemisphere(const Vec3& normal, Sampler* sampler)
	{
		// Sample uniformly on a unit hemisphere
		float u1 = sampler->next();
		float u2 = sampler->next();

		float r = sqrt(1.0f - u1 * u1);
		float phi = 2.0f * 3.14159f * u2;

		Vec3 localDir(r * cos(phi), r * sin(phi), u1);

		// Create a coordinate system based on the normal
		Vec3 tangent, bitangent;

		// Find a non-parallel vector to create a coordinate system
		if (fabs(normal.x) > fabs(normal.y))
			tangent = Vec3(normal.z, 0, -normal.x).normalize();
		else
			tangent = Vec3(0, -normal.z, normal.y).normalize();

		bitangent = Cross(normal, tangent);

		// Transform to world space
		return tangent * localDir.x + bitangent * localDir.y + normal * localDir.z;
	}

	void lightTracePath(Ray& r, Colour pathThroughput, Colour Le, Sampler* sampler)
	{
		// Similar to path tracing but no direct lighting
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);

		if (shadingData.t < FLT_MAX)
		{
			// Stop if we hit a light source
			if (shadingData.bsdf->isLight())
				return;

			// Connect this intersection point to the camera
			connectToCamera(shadingData.x, shadingData.sNormal, pathThroughput);

			// Russian roulette for path termination
			float russianRouletteProbability = min(pathThroughput.Lum(), 0.9f);
			if (sampler->next() >= russianRouletteProbability)
				return;

			pathThroughput = pathThroughput / russianRouletteProbability;

			// Sample BSDF to continue the path
			Colour bsdfValue;
			float pdf;
			Vec3 wi = shadingData.bsdf->sample(shadingData, sampler, bsdfValue, pdf);

			if (pdf <= 0.0f)
				return;

			pathThroughput = pathThroughput * bsdfValue * fabsf(Dot(wi, shadingData.sNormal)) / pdf;

			// Create the next ray
			Ray nextRay;
			nextRay.init(shadingData.x + (wi * EPSILON), wi);

			// Continue the light path recursively
			lightTracePath(nextRay, pathThroughput, Le, sampler);
		}
	}



	void render()
	{
		film->incrementSPP();

		int rowsPerThread = film->height / numProcs;

		// First pass: standard rendering and collect data for denoising
		auto renderChunk = [this](int startY, int endY, int threadIdx) {
			MTRandom* sampler = &samplers[threadIdx];

			// Mix of path tracing and light tracing
			// Using 80% path tracing, 20% light tracing as an example balance
			float lightTraceRatio = 0.0f; // Adjust as needed

			// Do some light tracing for this thread
			int lightTraceSamples = (int)(lightTraceRatio * (endY - startY) * film->width);
			for (int i = 0; i < lightTraceSamples; i++) {
				lightTrace(sampler);
			}

			// Now do path tracing for pixels
			for (unsigned int y = startY; y < endY; y++)
			{
				for (unsigned int x = 0; x < film->width; x++)
				{
					// Only do path tracing for the remaining percentage of samples
					if (sampler->next() >= lightTraceRatio) {
						float px = x + 0.5f;
						float py = y + 0.5f;
						Ray ray = scene->camera.generateRay(px, py);

						// Get color
						Colour pathThroughput = Colour(1.0f, 1.0f, 1.0f);
						int depth = 0;
						Colour col = pathTrace(ray, pathThroughput, depth, sampler, true);

						// Store values for denoiser
						int pixelIndex = (y * film->width + x) * 3;

						// Store color for denoiser
						colorBuffer[pixelIndex] = col.r;
						colorBuffer[pixelIndex + 1] = col.g;
						colorBuffer[pixelIndex + 2] = col.b;

						// Get and store albedo
						Colour albedoCol = albedo(ray);
						albedoBuffer[pixelIndex] = albedoCol.r;
						albedoBuffer[pixelIndex + 1] = albedoCol.g;
						albedoBuffer[pixelIndex + 2] = albedoCol.b;

						// Get and store normals
						Colour normalCol = viewNormals(ray);
						normalBuffer[pixelIndex] = normalCol.r;
						normalBuffer[pixelIndex + 1] = normalCol.g;
						normalBuffer[pixelIndex + 2] = normalCol.b;

						// Normal rendering path (still needed for film accumulation)
						film->splat(px, py, col);
						unsigned char r = (unsigned char)(col.r * 255);
						unsigned char g = (unsigned char)(col.g * 255);
						unsigned char b = (unsigned char)(col.b * 255);
						film->tonemap(x, y, r, g, b);
						canvas->draw(x, y, r, g, b);
					}
				}
			}
			};

		for (int i = 0; i < numProcs; i++)
		{
			int startY = i * rowsPerThread;
			int endY = (i == numProcs - 1) ? film->height : (i + 1) * rowsPerThread;
			threads[i] = new std::thread(renderChunk, startY, endY, i);
		}

		for (int i = 0; i < numProcs; i++)
		{
			threads[i]->join();
			delete threads[i];
			threads[i] = nullptr;
		}

		// Apply denoising if enabled
		if (enableDenoising) {
			applyDenoising();
		}
	}

	// Run the denoiser and update the display
	void applyDenoising()
	{
		// Only run denoiser if strength > 0
		if (denoiseStrength > 0.0f) {
			// Step 1: Run the denoiser on the collected buffers
			denoise(colorBuffer, albedoBuffer, normalBuffer, outputBuffer, film->width, film->height);

			// Step 2: Copy the film buffer to preserve accumulated samples
			Colour* tempFilm = new Colour[film->width * film->height];
			memcpy(tempFilm, film->film, film->width * film->height * sizeof(Colour));

			// Step 3: Clear the film for reuse
			Film* tempFilmObject = film;
			film = new Film();
			film->init((unsigned int)scene->camera.width, (unsigned int)scene->camera.height, tempFilmObject->filter);
			film->SPP = tempFilmObject->SPP;  // Preserve the SPP count

			// Step 4: Apply the blended result to the screen and resplat to film
			for (unsigned int y = 0; y < film->height; y++) {
				for (unsigned int x = 0; x < film->width; x++) {
					int pixelIndex = (y * film->width + x) * 3;
					int filmIndex = y * film->width + x;

					// Get the denoised pixel
					float r_denoised = outputBuffer[pixelIndex];
					float g_denoised = outputBuffer[pixelIndex + 1];
					float b_denoised = outputBuffer[pixelIndex + 2];

					// Get the original pixel
					float r_original = colorBuffer[pixelIndex];
					float g_original = colorBuffer[pixelIndex + 1];
					float b_original = colorBuffer[pixelIndex + 2];

					// Blend based on denoiseStrength
					float r = r_original * (1.0f - denoiseStrength) + r_denoised * denoiseStrength;
					float g = g_original * (1.0f - denoiseStrength) + g_denoised * denoiseStrength;
					float b = b_original * (1.0f - denoiseStrength) + b_denoised * denoiseStrength;

					// Create a color object
					Colour blendedCol(r, g, b);

					// Splat it to the film (at exact pixel centers to avoid filtering)
					float px = x + 0.5f;
					float py = y + 0.5f;
					film->splat(px, py, blendedCol * (float)film->SPP);

					// Update the canvas
					unsigned char rByte, gByte, bByte;
					film->tonemap(x, y, rByte, gByte, bByte);
					canvas->draw(x, y, rByte, gByte, bByte);
				}
			}

			// Clean up temporary objects
			delete tempFilmObject;
			delete[] tempFilm;
		}
		else {
			// If denoiseStrength is 0, just display the original image
			// (The original image is already displayed by the render function)
		}
	}

	// OIDN denoising function
	void denoise(float* colorBuffer, float* albedoBuffer, float* normalBuffer, float* outputBuffer, int width, int height)
	{
		// Create device
		oidn::DeviceRef device = oidn::newDevice();
		device.commit();

		// Create properly aligned buffers on the device
		oidn::BufferRef colorBufferOIDN = device.newBuffer(width * height * 3 * sizeof(float));
		oidn::BufferRef albedoBufferOIDN = device.newBuffer(width * height * 3 * sizeof(float));
		oidn::BufferRef normalBufferOIDN = device.newBuffer(width * height * 3 * sizeof(float));
		oidn::BufferRef outputBufferOIDN = device.newBuffer(width * height * 3 * sizeof(float));

		// Copy data to device buffers
		memcpy(colorBufferOIDN.getData(), colorBuffer, width * height * 3 * sizeof(float));
		memcpy(albedoBufferOIDN.getData(), albedoBuffer, width * height * 3 * sizeof(float));
		memcpy(normalBufferOIDN.getData(), normalBuffer, width * height * 3 * sizeof(float));

		// Create and configure the filter
		oidn::FilterRef filter = device.newFilter("RT");
		filter.setImage("color", colorBufferOIDN, oidn::Format::Float3, width, height);
		filter.setImage("albedo", albedoBufferOIDN, oidn::Format::Float3, width, height);
		filter.setImage("normal", normalBufferOIDN, oidn::Format::Float3, width, height);
		filter.setImage("output", outputBufferOIDN, oidn::Format::Float3, width, height);
		filter.set("hdr", true);
		filter.commit();

		// Execute the filter
		filter.execute();

		// Check for errors
		const char* errorMessage;
		if (device.getError(errorMessage) != oidn::Error::None) {
			printf("OIDN Error: %s\n", errorMessage);
		}

		// Copy the result back to the output buffer
		memcpy(outputBuffer, outputBufferOIDN.getData(), width * height * 3 * sizeof(float));
	}

	// Toggle denoising on/off
	void toggleDenoising(bool enable) {
		enableDenoising = enable;
	}

	// Destructor for clean up
	~RayTracer() {
		delete[] samplers;
		delete[] threads;
		delete film;

		// Clean up denoiser buffers
		delete[] colorBuffer;
		delete[] albedoBuffer;
		delete[] normalBuffer;
		delete[] outputBuffer;
	}

	int getSPP()
	{
		return film->SPP;
	}

	void saveHDR(std::string filename)
	{
		film->save(filename);
	}

	void savePNG(std::string filename)
	{
		stbi_write_png(filename.c_str(), canvas->getWidth(), canvas->getHeight(), 3, canvas->getBackBuffer(), canvas->getWidth() * 3);
	}
};