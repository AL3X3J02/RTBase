#pragma once

#include<algorithm>
#include "Core.h"
#include "Sampling.h"
#define EPSILON 0.000001f

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

class Triangle
{
public:
	Vertex vertices[3];
	Vec3 e1; // Edge 1
	Vec3 e2; // Edge 2
	Vec3 n; // Geometric Normal
	float area; // Triangle area
	float invarea; // Inverse area
	float d; // For ray triangle if needed
	unsigned int materialIndex;
	void init(Vertex v0, Vertex v1, Vertex v2, unsigned int _materialIndex)
	{
		materialIndex = _materialIndex;
		vertices[0] = v0;
		vertices[1] = v1;
		vertices[2] = v2;
		e1 = vertices[1].p - vertices[0].p;
		e2 = vertices[2].p - vertices[0].p;
		n = e1.cross(e2).normalize();
		area = e1.cross(e2).length() * 0.5f;
		invarea = 1.0f / area;
		d = -Dot(n, vertices[0].p);
	}
	Vec3 centre() const
	{
		return (vertices[0].p + vertices[1].p + vertices[2].p) / 3.0f;
	}

	// Möller-Trumbore
	bool rayIntersect(const Ray& r, float& t, float& u, float& v) const
	{
		Vec3 s = r.o - vertices[0].p;
		Vec3 s1 = r.dir.cross(e2);
		Vec3 s2 = s.cross(e1);
		float invDet = s1.dot(e1);
		if (std::abs(invDet) < EPSILON) return false;
		invDet = 1.0f / invDet;

		t = s2.dot(e2) * invDet;
		u = s1.dot(s) * invDet;
		v = s2.dot(r.dir) * invDet;
		if (t < 0.0f || u < 0.0f || v < 0.0f || u + v > 1.0f) return false;

		return true;

	}
	void interpolateAttributes(const float alpha, const float beta, const float gamma, Vec3& interpolatedNormal, float& interpolatedU, float& interpolatedV) const
	{
		interpolatedNormal = vertices[0].normal * alpha + vertices[1].normal * beta + vertices[2].normal * gamma;
		interpolatedNormal = interpolatedNormal.normalize();
		interpolatedU = vertices[0].u * alpha + vertices[1].u * beta + vertices[2].u * gamma;
		interpolatedV = vertices[0].v * alpha + vertices[1].v * beta + vertices[2].v * gamma;
	}

	Vec3 sample(Sampler* sampler, float& pdf)
	{
		float u1 = sampler->next();
		float u2 = sampler->next();

		float su1 = sqrt(u1);
		float alpha = 1.0f - su1;
		float beta = u2 * su1;
		float gamma = 1.0f - alpha - beta;

		Vec3 p = vertices[0].p * alpha + vertices[1].p * beta + vertices[2].p * gamma;
		pdf = invarea;

		return p;
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
	unsigned int offset = 0;
	unsigned int num = 0;

	bool isLeaf()
	{
		return (r == NULL && l == NULL);
	}

	BVHNode()
	{
		r = NULL;
		l = NULL;
	}
	~BVHNode()
	{
		if (r) delete r;
		if (l) delete l;
	}
	// Note there are several options for how to implement the build method. Update this as required
	void build(std::vector<Triangle>& inputTriangles, int start, int end)
	{

		// Add BVH building code here
		// Calculate bounds
		for (int i = start; i < end; i++) {
			bounds.extend(inputTriangles[i].vertices[0].p);
			bounds.extend(inputTriangles[i].vertices[1].p);
			bounds.extend(inputTriangles[i].vertices[2].p);
		}

		// If it has less than 8 triangles, it is a leaf node
		int numTriangles = end - start;
		if (numTriangles <= MAXNODE_TRIANGLES) {
			offset = start;
			num = numTriangles;
			return;
		}

		//SAH
		float bestCost = FLT_MAX;
		int bestAxis = -1;
		int bestSplit = -1;
		//std::vector<Triangle> originalTriangles(inputTriangles.begin() + start, inputTriangles.begin() + end);
		for (int axis = 0; axis < 3; axis++) {
			//std::vector<Triangle> iTriangles(originalTriangles);
			std::vector<Triangle> tempTriangles(inputTriangles.begin() + start, inputTriangles.begin() + end);

			std::stable_sort(tempTriangles.begin(), tempTriangles.end(),
				[axis](const Triangle& a, const Triangle& b) {
					float aCenter = a.vertices[0].p.coords[axis] + a.vertices[1].p.coords[axis] + a.vertices[2].p.coords[axis];
					float bCenter = b.vertices[0].p.coords[axis] + b.vertices[1].p.coords[axis] + b.vertices[2].p.coords[axis];
					return aCenter < bCenter;
				});

			std::vector<AABB> leftBounds(numTriangles), rightBounds(numTriangles);
			AABB box1, box2;
			for (int i = 0; i < numTriangles; i++) {
				box1.extend(tempTriangles[i].vertices[0].p);
				box1.extend(tempTriangles[i].vertices[1].p);
				box1.extend(tempTriangles[i].vertices[2].p);
				leftBounds[i] = box1;
			}
			for (int i = numTriangles - 1; i >= 0; i--) {
				box2.extend(tempTriangles[i].vertices[0].p);
				box2.extend(tempTriangles[i].vertices[1].p);
				box2.extend(tempTriangles[i].vertices[2].p);
				rightBounds[i] = box2;
			}

			// Calculate the best split
			for (int i = 1; i < numTriangles; i++) {

				float cost = TRAVERSE_COST +
					leftBounds[i - 1].area() / bounds.area() * i * TRIANGLE_COST +
					rightBounds[i].area() / bounds.area() * (numTriangles - i) * TRIANGLE_COST;

				if (cost < bestCost) {
					bestCost = cost;
					bestAxis = axis;
					bestSplit = start + i;
				}
			}
		}
		// can not find best spilt
		if (bestAxis == -1 || bestSplit <= start || bestSplit >= end)
		{
			offset = start;
			num = numTriangles;
			return;
		}

		//sort triangle vector according to best axis
		std::stable_sort(inputTriangles.begin() + start, inputTriangles.begin() + end,
			[bestAxis](const Triangle& a, const Triangle& b) {
				float aCenter = a.vertices[0].p.coords[bestAxis] + a.vertices[1].p.coords[bestAxis] + a.vertices[2].p.coords[bestAxis];
				float bCenter = b.vertices[0].p.coords[bestAxis] + b.vertices[1].p.coords[bestAxis] + b.vertices[2].p.coords[bestAxis];
				return aCenter < bCenter;
			});

		l = new BVHNode();
		r = new BVHNode();
		l->build(inputTriangles, start, bestSplit);
		r->build(inputTriangles, bestSplit, end);
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
		if (isLeaf())
		{
			// Test only triangles in this leaf node's range
			for (int i = offset; i < offset + num; i++)
			{
				float triT, u, v;
				if (triangles[i].rayIntersect(ray, triT, u, v))
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
		if (isLeaf())
		{
			// Test only triangles in this leaf node's range
			for (int i = offset; i < offset + num; i++)
			{
				const Triangle& tri = triangles[i];
				float triT, u, v;

				// Use the Möller-Trumbore algorithm for ray-triangle intersection
				if (tri.rayIntersect(ray, triT, u, v))
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



