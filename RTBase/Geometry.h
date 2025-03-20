#pragma once
#define EPSILON 0.0001f
#include "Core.h"
#include "Sampling.h"

class Ray
{
public:
	Vec3 o;
	Vec3 dir;
	Vec3 invDir;
	Ray()
	{
	}
	Ray(Vec3 _o, Vec3 _d)
	{
		init(_o, _d);
	}
	void init(Vec3 _o, Vec3 _d)
	{
		o = _o;
		dir = _d;
		invDir = Vec3(1.0f / dir.x, 1.0f / dir.y, 1.0f / dir.z);
	}
	Vec3 at(const float t) const
	{
		return (o + (dir * t));
	}
};

class Plane
{
public:
	Vec3 n;
	float d;

	void init(const Vec3& _n, float _d)
	{
		n = _n.normalize();
		d = _d;
	}

	void initFromPoints(const Vec3& p0, const Vec3& p1, const Vec3& p2)
	{
		Vec3 e1 = p1 - p0;
		Vec3 e2 = p2 - p0;
		n = Cross(e1, e2).normalize();
		d = Dot(n, p0);
	}

	// Ray-plane intersection
	bool rayIntersect(Ray& r, float& t)
	{
		float denom = Dot(n, r.dir);

		if (std::abs(denom) < EPSILON)
			return false;

		t = (d - Dot(n, r.o)) / denom;

		return t >= 0.f;
	}
};

#define EPSILON 0.001f

const static float g_fEPSILON = 0.000001f;
class Triangle
{
public:
	Vertex vertices[3];
	Vec3 e1; // Edge 1
	Vec3 e2; // Edge 2
	Vec3 n; // Geometric Normal

	float area; // Triangle area
	float d; // For ray triangle if needed
	unsigned int materialIndex;

	Vec3 v0v1;
	Vec3 v0v2;

	void init(const Vertex& v0, const Vertex& v1, const Vertex& v2, const unsigned int& _materialIndex)
	{
		materialIndex = _materialIndex;
		vertices[0] = v0;
		vertices[1] = v1;
		vertices[2] = v2;

		e1 = vertices[2].p - vertices[1].p;
		e2 = vertices[0].p - vertices[2].p;
		n = e1.cross(e2).normalize();
		area = e1.cross(e2).length() * 0.5f;
		d = Dot(n, vertices[0].p);

		v0v1 = v1.p - v0.p;
		v0v2 = v2.p - v0.p;
	}
	Vec3 centre() const
	{
		return (vertices[0].p + vertices[1].p + vertices[2].p) / 3.0f;
	}
	// Add code here
	bool rayIntersect(const Ray& r, float& t, float& u, float& v) const
	{
		float denom = Dot(n, r.dir);
		if (denom == 0) { return false; }
		t = (d - Dot(n, r.o)) / denom;
		if (t < 0) { return false; }
		Vec3 p = r.at(t);
		float invArea = 1.0f / Dot(e1.cross(e2), n);
		u = Dot(e1.cross(p - vertices[1].p), n) * invArea;
		if (u < 0 || u > 1.0f) { return false; }
		v = Dot(e2.cross(p - vertices[2].p), n) * invArea;
		if (v < 0 || (u + v) > 1.0f) { return false; }
		return true;
	}
	//Add Moller trumbore
	bool rayIntersectMollerTrumbore(const Ray& r, float& t, float& u, float& v) const
	{
		const Vec3 pvec = Cross(r.dir, v0v2);
		const float det = Dot(v0v1, pvec);

		if (det < g_fEPSILON
			&& det > -g_fEPSILON)
			return false;

		const float invDet = 1.0f / det;
		const Vec3 tvec = r.o - vertices[0].p;
		u = Dot(tvec, pvec) * invDet;

		if (u < 0.0f || u > 1.0f)
			return false;

		const Vec3 qvec = Cross(tvec, v0v1);
		v = Dot(r.dir, qvec) * invDet;

		if (v < 0.0f || (u + v) > 1.0f)
			return false;

		t = Dot(v0v2, qvec) * invDet;

		return true;
	}
	/*bool rayIntersect(const Ray& r, float& t, float& u, float& v) const
	{
		t = (d-Dot(n, r.o) ) / Dot(n, r.dir);
		if (t < 0.f) return false;

		Vec3 P = r.at(t);
		Vec3 q1 = P - vertices[0].p;
		Vec3 C1 = Cross(e1, q1);

		u = Dot(C1 , n) / area;
		if (u > 1.f || u < 0.f)return false;

		Vec3 q2 = P - vertices[1].p;
		Vec3 C2 = Cross(e2, q2);

		v = Dot(C2, n) / area;
		if (v > 1.f || v < 0.f) return false;

		if (u + v > 1.f) return false;

		return true;
	}*/
	void interpolateAttributes(const float alpha, const float beta, const float gamma, Vec3& interpolatedNormal, float& interpolatedU, float& interpolatedV) const
	{
		interpolatedNormal = vertices[0].normal * alpha + vertices[1].normal * beta + vertices[2].normal * gamma;
		interpolatedNormal = interpolatedNormal.normalize();
		interpolatedU = vertices[0].u * alpha + vertices[1].u * beta + vertices[2].u * gamma;
		interpolatedV = vertices[0].v * alpha + vertices[1].v * beta + vertices[2].v * gamma;
	}
	// Add code here

	Vec3 sample(Sampler* sampler, float& pdf)
	{
		float r1 = sampler->next();
		float r2 = sampler->next();

		float alpha = 1.0f - sqrt(r1);
		float beta = r2 * sqrt(r1);
		float gamma = 1.0f - alpha - beta;

		Vec3 point = vertices[0].p * alpha + vertices[1].p * beta + vertices[2].p * gamma;

		pdf = 1.0f / area;

		return point;
	}
	Vec3 gNormal()
	{
		return (n * (Dot(vertices[0].normal, n) > 0 ? 1.0f : -1.0f));
	}
};

class AABB
{
public:
	Vec3 max;
	Vec3 min;

	AABB()
	{
		reset();
	}

	void reset()
	{
		max = Vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
		min = Vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	}

	void extend(const Vec3 p)
	{
		max = Vec3(std::max(max.x, p.x), std::max(max.y, p.y), std::max(max.z, p.z));
		min = Vec3(std::min(min.x, p.x), std::min(min.y, p.y), std::min(min.z, p.z));
	}

	// Ray-AABB intersection without returning t
	bool rayAABB(const Ray& r)
	{
		float tmin = 0.0f;
		float tmax = FLT_MAX;

		float invDirX = r.invDir.x;
		float t1x = (min.x - r.o.x) * invDirX;
		float t2x = (max.x - r.o.x) * invDirX;
		if (invDirX < 0.0f) std::swap(t1x, t2x);
		tmin = std::max(t1x, tmin);
		tmax = std::min(t2x, tmax);
		if (tmin > tmax) return false;

		float invDirY = r.invDir.y;
		float t1y = (min.y - r.o.y) * invDirY;
		float t2y = (max.y - r.o.y) * invDirY;
		if (invDirY < 0.0f) std::swap(t1y, t2y);
		tmin = std::max(t1y, tmin);
		tmax = std::min(t2y, tmax);
		if (tmin > tmax) return false;

		float invDirZ = r.invDir.z;
		float t1z = (min.z - r.o.z) * invDirZ;
		float t2z = (max.z - r.o.z) * invDirZ;
		if (invDirZ < 0.0f) std::swap(t1z, t2z);
		tmin = std::max(t1z, tmin);
		tmax = std::min(t2z, tmax);
		if (tmin > tmax) return false;

		return tmin >= 0.0f;
	}

	// Ray-AABB intersection with returning t
	bool rayAABB(const Ray& r, float& t)
	{
		float tmin = 0.0f;
		float tmax = FLT_MAX;

		float invDirX = r.invDir.x;
		float t1x = (min.x - r.o.x) * invDirX;
		float t2x = (max.x - r.o.x) * invDirX;
		if (invDirX < 0.0f) std::swap(t1x, t2x);
		tmin = std::max(t1x, tmin);
		tmax = std::min(t2x, tmax);
		if (tmin > tmax) return false;

		float invDirY = r.invDir.y;
		float t1y = (min.y - r.o.y) * invDirY;
		float t2y = (max.y - r.o.y) * invDirY;
		if (invDirY < 0.0f) std::swap(t1y, t2y);
		tmin = std::max(t1y, tmin);
		tmax = std::min(t2y, tmax);
		if (tmin > tmax) return false;

		float invDirZ = r.invDir.z;
		float t1z = (min.z - r.o.z) * invDirZ;
		float t2z = (max.z - r.o.z) * invDirZ;
		if (invDirZ < 0.0f) std::swap(t1z, t2z);
		tmin = std::max(t1z, tmin);
		tmax = std::min(t2z, tmax);
		if (tmin > tmax) return false;

		if (tmin >= 0.0f)
		{
			t = tmin;
			return true;
		}

		return false;
	}

	float area()
	{
		Vec3 size = max - min;
		return 2.0f * (size.x * size.y + size.y * size.z + size.x * size.z);
	}
};

class Sphere
{
public:
	Vec3 centre;
	float radius;

	void init(const Vec3& _centre, float _radius)
	{
		centre = _centre;
		radius = _radius;
	}

	// Ray-sphere intersection
	bool rayIntersect(Ray& r, float& t)
	{
		Vec3 oc = r.o - centre;
		float a = Dot(r.dir, r.dir);
		float b = 2.0f * Dot(oc, r.dir);
		float c = Dot(oc, oc) - radius * radius;

		float discriminant = b * b - 4 * a * c;

		if (discriminant < 0.0f)
			return false;

		float sqrtDiscriminant = std::sqrt(discriminant);
		float t0 = (-b - sqrtDiscriminant) / (2.0f * a);
		float t1 = (-b + sqrtDiscriminant) / (2.0f * a);

		if (t0 > 0.0f && t1 > 0.0f)
			t = std::min(t0, t1);
		else if (t0 > 0.0f)
			t = t0;
		else if (t1 > 0.0f)
			t = t1;
		else
			return false;

		return true;
	}
};

struct IntersectionData
{
	unsigned int ID;
	float t;
	float alpha;
	float beta;
	float gamma;
};

#define MAXNODE_TRIANGLES 8
#define TRAVERSE_COST 1.0f
#define TRIANGLE_COST 2.0f
#define BUILD_BINS 32

class BVHNode
{
public:
	AABB bounds;
	BVHNode* r;
	BVHNode* l;
	// This can store an offset and number of triangles in a global triangle list for example
	// But you can store this however you want!
	// unsigned int offset;
	// unsigned char num;
	int startIndex;
	int triangleCount;

	BVHNode()
	{
		r = NULL;
		l = NULL;
	}
	// Note there are several options for how to implement the build method. Update this as required
	void build(std::vector<Triangle>& inputTriangles, std::vector<Triangle>& outputTriangles, int depth = 0)
	{
		// Record the starting index in the output array
		startIndex = outputTriangles.size();

		// If the input is empty, return without doing anything
		if (inputTriangles.empty())
			return;

		// Calculate bounds for all triangles
		bounds.reset();
		for (const Triangle& tri : inputTriangles)
		{
			bounds.extend(tri.vertices[0].p);
			bounds.extend(tri.vertices[1].p);
			bounds.extend(tri.vertices[2].p);
		}

		// If we have a small number of triangles, make this a leaf node
		if (inputTriangles.size() <= MAXNODE_TRIANGLES || depth > 20) // Add depth limit
		{
			triangleCount = inputTriangles.size();
			for (const Triangle& tri : inputTriangles)
				outputTriangles.push_back(tri);
			return;
		}

		// Find the axis with the largest extent
		Vec3 extents = bounds.max - bounds.min;
		int axis = 0;
		if (extents.y > extents.x) axis = 1;
		if (extents.z > (axis == 0 ? extents.x : extents.y)) axis = 2;

		// SAH evaluation
		float bestCost = FLT_MAX;
		int bestBin = -1;

		// Calculate min and max for the axis we're splitting on
		float axisMin = (axis == 0) ? bounds.min.x : ((axis == 1) ? bounds.min.y : bounds.min.z);
		float axisMax = (axis == 0) ? bounds.max.x : ((axis == 1) ? bounds.max.y : bounds.max.z);
		float axisRange = axisMax - axisMin;

		struct Bin {
			AABB bounds;
			int count = 0;
		};
		std::vector<Bin> bins(BUILD_BINS);

		// Assign triangles to bins based on their centers
		for (const Triangle& tri : inputTriangles)
		{
			Vec3 center = tri.centre();
			float centroid;
			if (axis == 0) centroid = center.x;
			else if (axis == 1) centroid = center.y;
			else centroid = center.z;

			int binIndex = std::min(BUILD_BINS - 1, static_cast<int>((centroid - axisMin) / axisRange * BUILD_BINS));

			bins[binIndex].count++;
			bins[binIndex].bounds.extend(tri.vertices[0].p);
			bins[binIndex].bounds.extend(tri.vertices[1].p);
			bins[binIndex].bounds.extend(tri.vertices[2].p);
		}

		// Evaluate SAH for each possible split position
		std::vector<AABB> leftAcc(BUILD_BINS);
		std::vector<AABB> rightAcc(BUILD_BINS);
		std::vector<int> leftCount(BUILD_BINS, 0);
		std::vector<int> rightCount(BUILD_BINS, 0);

		// Left-to-right sweep to accumulate left bounds and counts
		leftAcc[0] = bins[0].bounds;
		leftCount[0] = bins[0].count;
		for (int i = 1; i < BUILD_BINS; i++)
		{
			leftAcc[i] = leftAcc[i - 1];
			if (bins[i].count > 0)
			{
				leftAcc[i].extend(bins[i].bounds.min);
				leftAcc[i].extend(bins[i].bounds.max);
			}
			leftCount[i] = leftCount[i - 1] + bins[i].count;
		}

		// Right-to-left sweep to accumulate right bounds and counts
		rightAcc[BUILD_BINS - 1] = bins[BUILD_BINS - 1].bounds;
		rightCount[BUILD_BINS - 1] = bins[BUILD_BINS - 1].count;
		for (int i = BUILD_BINS - 2; i >= 0; i--)
		{
			rightAcc[i] = rightAcc[i + 1];
			if (bins[i].count > 0)
			{
				rightAcc[i].extend(bins[i].bounds.min);
				rightAcc[i].extend(bins[i].bounds.max);
			}
			rightCount[i] = rightCount[i + 1] + bins[i].count;
		}

		// Evaluate SAH for each potential split
		for (int i = 0; i < BUILD_BINS - 1; i++)
		{
			if (leftCount[i] == 0 || rightCount[i + 1] == 0)
				continue;

			float surfaceAreaLeft = leftAcc[i].area();
			float surfaceAreaRight = rightAcc[i + 1].area();
			float totalSurfaceArea = bounds.area();

			float costLeft = (surfaceAreaLeft / totalSurfaceArea) * leftCount[i] * TRIANGLE_COST;
			float costRight = (surfaceAreaRight / totalSurfaceArea) * rightCount[i + 1] * TRIANGLE_COST;
			float cost = TRAVERSE_COST + costLeft + costRight;

			if (cost < bestCost)
			{
				bestCost = cost;
				bestBin = i;
			}
		}

		// Calculate the cost of not splitting (leaf cost)
		float leafCost = inputTriangles.size() * TRIANGLE_COST;

		if (bestCost >= leafCost || bestBin == -1)
		{
			for (const Triangle& tri : inputTriangles)
				outputTriangles.push_back(tri);
			return;
		}

		// Split position
		float splitPos = axisMin + (bestBin + 1) * axisRange / BUILD_BINS;

		// Partition triangles
		std::vector<Triangle> leftTriangles;
		std::vector<Triangle> rightTriangles;

		for (const Triangle& tri : inputTriangles)
		{
			Vec3 center = tri.centre();
			float value;
			if (axis == 0) value = center.x;
			else if (axis == 1) value = center.y;
			else value = center.z;

			if (value <= splitPos)
				leftTriangles.push_back(tri);
			else
				rightTriangles.push_back(tri);
		}

		// Handle edge cases
		if (leftTriangles.empty() || rightTriangles.empty())
		{
			size_t halfSize = inputTriangles.size() / 2;
			leftTriangles.clear();
			rightTriangles.clear();

			for (size_t i = 0; i < inputTriangles.size(); i++)
			{
				if (i < halfSize)
					leftTriangles.push_back(inputTriangles[i]);
				else
					rightTriangles.push_back(inputTriangles[i]);
			}
		}

		l = new BVHNode();
		r = new BVHNode();

		l->build(leftTriangles, outputTriangles, depth + 1);
		r->build(rightTriangles, outputTriangles, depth + 1);
	}
	void traverse(const Ray& ray, const std::vector<Triangle>& triangles, IntersectionData& intersection)
	{
		// First check if ray intersects bounding box
		float t;
		if (!bounds.rayAABB(ray, t))
			return;

		// If we've already found a closer intersection, skip this node
		if (t > intersection.t)
			return;

		// If this is a leaf node (no children)
		if (l == NULL && r == NULL)
		{
			// Test only triangles in this leaf node's range
			for (int i = startIndex; i < startIndex + triangleCount; i++)
			{
				float triT, u, v;
				if (triangles[i].rayIntersectMollerTrumbore(ray, triT, u, v))
				{
					if (triT < intersection.t)
					{
						intersection.t = triT;
						intersection.ID = i;
						intersection.alpha = u;
						intersection.beta = v;
						intersection.gamma = 1.0f - (u + v);
					}
				}
			}
		}
		else
		{
			// Not a leaf node, so traverse the children
			// Traverse the closer child first (front-to-back traversal)
			BVHNode* firstNode = l;
			BVHNode* secondNode = r;

			// Determine which child node is closer to the ray origin
			float tLeft = FLT_MAX;
			float tRight = FLT_MAX;

			if (l) l->bounds.rayAABB(ray, tLeft);
			if (r) r->bounds.rayAABB(ray, tRight);

			// Swap the traversal order if the right child is closer
			if (tRight < tLeft)
			{
				firstNode = r;
				secondNode = l;
			}

			// Traverse the closer child first
			if (firstNode) firstNode->traverse(ray, triangles, intersection);

			// Only traverse the second child if we could potentially find a closer intersection
			if (secondNode && intersection.t > t) secondNode->traverse(ray, triangles, intersection);
		}
	}
	IntersectionData traverse(const Ray& ray, const std::vector<Triangle>& triangles)
	{
		IntersectionData intersection;
		intersection.t = FLT_MAX;
		traverse(ray, triangles, intersection);
		return intersection;
	}
	bool traverseVisible(const Ray& ray, const std::vector<Triangle>& triangles, const float maxT)
	{
		// First check if the ray intersects the bounding box of this node
		float t;
		if (!bounds.rayAABB(ray, t))
			return true; // No intersection with bounds, the path is clear

		// If the intersection with bounds is beyond our max distance, the path is clear
		if (t > maxT)
			return true;

		// If this is a leaf node (no children)
		if (l == NULL && r == NULL)
		{
			// Test only triangles in this leaf node's range
			for (int i = startIndex; i < startIndex + triangleCount; i++)
			{
				const Triangle& tri = triangles[i];
				float triT, u, v;

				// Use the Möller-Trumbore algorithm for ray-triangle intersection
				if (tri.rayIntersectMollerTrumbore(ray, triT, u, v))
				{
					// If we hit something within our max distance, the path is blocked
					if (triT < maxT)
						return false;
				}
			}

			// No triangles were hit within our max distance
			return true;
		}
		else
		{
			// Not a leaf node, so traverse the children
			// Traverse the closer child first (front-to-back traversal)
			BVHNode* firstNode = l;
			BVHNode* secondNode = r;

			// Determine which child node is closer to the ray origin
			float tLeft = FLT_MAX;
			float tRight = FLT_MAX;

			if (l) l->bounds.rayAABB(ray, tLeft);
			if (r) r->bounds.rayAABB(ray, tRight);

			// Swap the traversal order if the right child is closer
			if (tRight < tLeft)
			{
				firstNode = r;
				secondNode = l;
			}

			// Traverse the closer child first
			if (firstNode && !firstNode->traverseVisible(ray, triangles, maxT))
				return false; // Path is blocked in the first child

			// Only traverse the second child if necessary
			if (secondNode && !secondNode->traverseVisible(ray, triangles, maxT))
				return false; // Path is blocked in the second child

			// No blocking geometry found in either child
			return true;
		}
	}
};



