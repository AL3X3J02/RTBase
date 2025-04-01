#pragma once

#include "Core.h"
#include "Geometry.h"
#include "Materials.h"
#include "Sampling.h"

#pragma warning( disable : 4244)

class SceneBounds
{
public:
    Vec3 sceneCentre;
    float sceneRadius;
};

class Light
{
public:
    virtual Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& emittedColour, float& pdf) = 0;
    virtual Colour evaluate(const ShadingData& shadingData, const Vec3& wi) = 0;
    virtual float PDF(const ShadingData& shadingData, const Vec3& wi) = 0;
    virtual bool isArea() = 0;
    virtual Vec3 normal(const ShadingData& shadingData, const Vec3& wi) = 0;
    virtual float totalIntegratedPower() = 0;
    virtual Vec3 samplePositionFromLight(Sampler* sampler, float& pdf) = 0;
    virtual Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf) = 0;
};

class AreaLight : public Light
{
public:
    Triangle* triangle = NULL;
    Colour emission;
    Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& emittedColour, float& pdf)
    {
        emittedColour = emission;
        return triangle->sample(sampler, pdf);
    }
    Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
    {
        if (Dot(wi, triangle->gNormal()) < 0)
        {
            return emission;
        }
        return Colour(0.0f, 0.0f, 0.0f);
    }
    float PDF(const ShadingData& shadingData, const Vec3& wi)
    {
        return 1.0f / triangle->area;
    }
    bool isArea()
    {
        return true;
    }
    Vec3 normal(const ShadingData& shadingData, const Vec3& wi)
    {
        return triangle->gNormal();
    }
    float totalIntegratedPower()
    {
        return (triangle->area * emission.Lum());
    }
    Vec3 samplePositionFromLight(Sampler* sampler, float& pdf)
    {
        return triangle->sample(sampler, pdf);
    }
    Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf)
    {
        // Add code to sample a direction from the light
        Vec3 wi = Vec3(0, 0, 1);
        pdf = 1.0f;
        Frame frame;
        frame.fromVector(triangle->gNormal());
        return frame.toWorld(wi);
    }
};

class BackgroundColour : public Light
{
public:
    Colour emission;
    BackgroundColour(Colour _emission)
    {
        emission = _emission;
    }
    Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
    {
        Vec3 wi = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
        pdf = SamplingDistributions::uniformSpherePDF(wi);
        reflectedColour = emission;
        return wi;
    }
    Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
    {
        return emission;
    }
    float PDF(const ShadingData& shadingData, const Vec3& wi)
    {
        return SamplingDistributions::uniformSpherePDF(wi);
    }
    bool isArea()
    {
        return false;
    }
    Vec3 normal(const ShadingData& shadingData, const Vec3& wi)
    {
        return -wi;
    }
    float totalIntegratedPower()
    {
        return emission.Lum() * 4.0f * M_PI;
    }
    Vec3 samplePositionFromLight(Sampler* sampler, float& pdf)
    {
        Vec3 p = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
        p = p * use<SceneBounds>().sceneRadius;
        p = p + use<SceneBounds>().sceneCentre;
        pdf = 4 * M_PI * use<SceneBounds>().sceneRadius * use<SceneBounds>().sceneRadius;
        return p;
    }
    Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf)
    {
        Vec3 wi = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
        pdf = SamplingDistributions::uniformSpherePDF(wi);
        return wi;
    }
};

// Slides 65-69 for light
// Fix for the EnvironmentMap class
// Fix for the EnvironmentMap class

// MIS HERE IMPLEMENT
class EnvironmentMap : public Light
{
public:
    Texture* env;
    // Grid for importance sampling
    float* cdf_v;               // Marginal CDF for v (rows)
    float** cdf_u;              // Conditional CDFs for u (columns) for each v
    float luminanceTotal;       // Total luminance for normalization
    bool samplingInitialized;   // Flag to check if sampling tables are initialized

    // Helper function to clamp values
    template <typename T>
    T clamp(T value, T min, T max) const
    {
        if (value < min) return min;
        if (value > max) return max;
        return value;
    }

    EnvironmentMap(Texture* _env)
    {
        env = _env;
        samplingInitialized = false;
        cdf_v = nullptr;
        cdf_u = nullptr;
        luminanceTotal = 0.0f;
        initImportanceSampling();
    }

    // Initialize the importance sampling grid
    void initImportanceSampling()
    {
        if (samplingInitialized || env == nullptr)
            return;

        int width = env->width;
        int height = env->height;

        // Allocate memory for distributions
        cdf_v = new float[height + 1];
        cdf_u = new float* [height];
        for (int i = 0; i < height; i++)
        {
            cdf_u[i] = new float[width + 1];
        }

        // Initialize with zeros
        cdf_v[0] = 0.0f;
        for (int i = 0; i < height; i++)
        {
            cdf_u[i][0] = 0.0f;
        }

        // Compute the marginal distributions p(v)
        float* row_luminance = new float[height];
        luminanceTotal = 0.0f;

        for (int i = 0; i < height; i++)
        {
            // Apply sinf weighting for spherical integration
            float sinTheta = sinf(((float)i + 0.5f) / (float)height * M_PI);
            row_luminance[i] = 0.0f;

            // Compute luminance for each row and build conditional CDFs
            for (int j = 0; j < width; j++)
            {
                float lum = env->texels[(i * width) + j].Lum() * sinTheta;
                row_luminance[i] += lum;
                cdf_u[i][j + 1] = cdf_u[i][j] + lum;
            }

            // Normalize the conditional CDF
            if (row_luminance[i] > 0.0f)
            {
                for (int j = 1; j <= width; j++)
                {
                    cdf_u[i][j] /= row_luminance[i];
                }
            }

            // Build the marginal CDF
            luminanceTotal += row_luminance[i];
            cdf_v[i + 1] = cdf_v[i] + row_luminance[i];
        }

        // Normalize the marginal CDF
        if (luminanceTotal > 0.0f)
        {
            for (int i = 1; i <= height; i++)
            {
                cdf_v[i] /= luminanceTotal;
            }
        }

        delete[] row_luminance;
        samplingInitialized = true;
    }

    // Binary search to find the index for a given target probability
    int findInterval(const float* cdf, int size, float target)
    {
        int low = 0, high = size - 1;
        while (low < high)
        {
            int mid = (low + high) >> 1;
            if (cdf[mid] <= target)
                low = mid + 1;
            else
                high = mid;
        }
        return clamp(low - 1, 0, size - 2);
    }

    Colour computeMIS(const ShadingData& shadingData, const Vec3& wi, bool isDirect = true)
    {
        if (!samplingInitialized)
        {
            // If sampling is not initialized, return the basic evaluation
            return evaluate(shadingData, wi);
        }

        // Convert direction to (u, v) texture coordinates
        float u = atan2f(wi.z, wi.x);
        u = (u < 0.0f) ? u + (2.0f * M_PI) : u;
        u = u / (2.0f * M_PI);
        float v = acosf(wi.y) / M_PI;

        // Convert to pixel coordinates
        int u_idx = clamp((int)(u * env->width), 0, env->width - 1);
        int v_idx = clamp((int)(v * env->height), 0, env->height - 1);

        // Number of samples (typically set to a fixed value or passed as a parameter)
        const int N = 1; // You might want to make this configurable

        // Calculate the weight based on the provided formulas
        float w;
        if (isDirect)
        {
            // Direct lighting MIS weight calculation
            float pA_x0 = 1.0f; // Probability of selecting the current point 
            float pA_x1_to_x0 = PDF(shadingData, wi); // Probability of sampling this direction
            w = pA_x0 / (pA_x0 + pA_x1_to_x0);
        }
        else
        {
            // Indirect lighting MIS weight calculation
            float pA_x1_to_x0 = 1.0f; // Probability of the previous path segment
            float pN_wi = PDF(shadingData, wi); // Probability of sampling this direction
            w = pA_x1_to_x0 / (pA_x1_to_x0 + pN_wi);
        }

        // Get the emission color
        Colour emissionColor = evaluate(shadingData, wi);

        // Apply MIS weight
        return emissionColor * w;
    }

    Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& emittedColour, float& pdf) override
    {
        if (!samplingInitialized)
        {
            initImportanceSampling();
        }

        int width = env->width;
        int height = env->height;

        // Sample a row using the marginal distribution
        float rv = sampler->next();
        int v_idx = findInterval(cdf_v, height + 1, rv);

        // Sample a column using the conditional distribution for the chosen row
        float ru = sampler->next();
        int u_idx = findInterval(cdf_u[v_idx], width + 1, ru);

        // Compute continuous sample coordinates (u, v) within selected cell
        float u_offset = (ru - cdf_u[v_idx][u_idx]) / (cdf_u[v_idx][u_idx + 1] - cdf_u[v_idx][u_idx]);
        float v_offset = (rv - cdf_v[v_idx]) / (cdf_v[v_idx + 1] - cdf_v[v_idx]);

        // Map to continuous coordinates in [0, 1)
        float u = (u_idx + u_offset) / width;
        float v = (v_idx + v_offset) / height;

        // Convert to spherical coordinates
        float theta = v * M_PI;
        float phi = u * 2.0f * M_PI;

        // Convert to direction vector
        Vec3 wi = sampleDirectionFromLight(sampler, pdf);

        // Use MIS for direct lighting calculation
        emittedColour = computeMIS(shadingData, wi, true);

        // Calculate the PDF
        pdf = PDF(shadingData, wi);

        return wi;
    }

    Vec3 sampleIndirect(const ShadingData& shadingData, Sampler* sampler, Colour& emittedColour, float& pdf)
{
    Vec3 wi = sampleDirectionFromLight(sampler, pdf);
    
    // Use MIS for indirect lighting calculation
    emittedColour = computeMIS(shadingData, wi, false);
    
    return wi;
}

    Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
    {
        float u = atan2f(wi.z, wi.x);
        u = (u < 0.0f) ? u + (2.0f * M_PI) : u;
        u = u / (2.0f * M_PI);
        float v = acosf(wi.y) / M_PI;
        return env->sample(u, v);
    }

    float PDF(const ShadingData& shadingData, const Vec3& wi)
    {
        if (!samplingInitialized)
            return SamplingDistributions::uniformSpherePDF(wi);

        // Convert direction to (u, v) texture coordinates
        float u = atan2f(wi.z, wi.x);
        u = (u < 0.0f) ? u + (2.0f * M_PI) : u;
        u = u / (2.0f * M_PI);
        float v = acosf(wi.y) / M_PI;

        // Convert to pixel coordinates
        int u_idx = clamp((int)(u * env->width), 0, env->width - 1);
        int v_idx = clamp((int)(v * env->height), 0, env->height - 1);

        // Get pixel luminance
        float sinTheta = sinf(v * M_PI);
        if (sinTheta == 0.0f)
            return 0.0f;

        float pixel_luminance = env->texels[(v_idx * env->width) + u_idx].Lum();

        // PDF = (pixel_luminance * sin(theta)) / (luminanceTotal * (2π * π) / (width * height))
        return (pixel_luminance * env->width * env->height) / (luminanceTotal * 2.0f * M_PI * M_PI * sinTheta);
    }

    bool isArea()
    {
        return false;
    }

    Vec3 normal(const ShadingData& shadingData, const Vec3& wi)
    {
        return -wi;
    }

    float totalIntegratedPower()
    {
        if (!samplingInitialized)
        {
            float total = 0;
            for (int i = 0; i < env->height; i++)
            {
                float st = sinf(((float)i / (float)env->height) * M_PI);
                for (int n = 0; n < env->width; n++)
                {
                    total += (env->texels[(i * env->width) + n].Lum() * st);
                }
            }
            total = total / (float)(env->width * env->height);
            return total * 4.0f * M_PI;
        }
        return luminanceTotal * 4.0f * M_PI;
    }

    // Move these functions outside of totalIntegratedPower
    Vec3 samplePositionFromLight(Sampler* sampler, float& pdf) override
    {
        // Samples a point on the bounding sphere of the scene
        Vec3 p = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
        p = p * use<SceneBounds>().sceneRadius;
        p = p + use<SceneBounds>().sceneCentre;
        pdf = 1.0f / (4 * M_PI * SQ(use<SceneBounds>().sceneRadius));
        return p;
    }

    Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf) override
    {
        // Using importance sampling if initialized, otherwise uniform sampling
        if (!samplingInitialized)
        {
            Vec3 wi = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
            pdf = SamplingDistributions::uniformSpherePDF(wi);
            return wi;
        }

        // Sample direction using the importance sampling tables
        int width = env->width;
        int height = env->height;

        float rv = sampler->next();
        int v_idx = findInterval(cdf_v, height + 1, rv);

        float ru = sampler->next();
        int u_idx = findInterval(cdf_u[v_idx], width + 1, ru);

        float u_offset = (ru - cdf_u[v_idx][u_idx]) / (cdf_u[v_idx][u_idx + 1] - cdf_u[v_idx][u_idx]);
        float v_offset = (rv - cdf_v[v_idx]) / (cdf_v[v_idx + 1] - cdf_v[v_idx]);

        float u = (u_idx + u_offset) / width;
        float v = (v_idx + v_offset) / height;

        float theta = v * M_PI;
        float phi = u * 2.0f * M_PI;

        Vec3 wi = Vec3(sinf(theta) * cosf(phi), cosf(theta), sinf(theta) * sinf(phi));

        float sinTheta = sinf(theta);
        if (sinTheta == 0.0f)
        {
            pdf = 0.0f;
        }
        else
        {
            float pixel_luminance = env->texels[(v_idx * width) + u_idx].Lum();
            pdf = (pixel_luminance * width * height) / (luminanceTotal * 2.0f * M_PI * M_PI * sinTheta);
        }

        return wi;
    }

    ~EnvironmentMap()
    {
        if (samplingInitialized && cdf_v != nullptr)
        {
            delete[] cdf_v;
            for (int i = 0; i < env->height; i++)
            {
                delete[] cdf_u[i];
            }
            delete[] cdf_u;
        }
    }
};