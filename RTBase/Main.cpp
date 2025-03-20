

#include "GEMLoader.h"
#include "Renderer.h"
#include "SceneLoader.h"
#define NOMINMAX
#include "GamesEngineeringBase.h"
#include <unordered_map>
#include <iostream>

void runTests()
{
    std::cout << "Running tests..." << std::endl;

    // Test AABB ray intersection
    {
        AABB aabb;
        aabb.min = Vec3(0, 0, 0);
        aabb.max = Vec3(1, 1, 1);

        // Test ray that intersects the AABB
        Ray ray1(Vec3(-1, 0.5, 0.5), Vec3(1, 0, 0)); // Ray pointing towards the AABB
        float t1;
        bool intersect1 = aabb.rayAABB(ray1, t1);
        std::cout << "Test 1 (Ray intersects AABB): " << (intersect1 ? "PASSED" : "FAILED") << std::endl;

        // Test ray that misses the AABB
        Ray ray2(Vec3(-1, 2, 2), Vec3(1, 0, 0)); // Ray pointing above the AABB
        float t2;
        bool intersect2 = aabb.rayAABB(ray2, t2);
        std::cout << "Test 2 (Ray misses AABB): " << (!intersect2 ? "PASSED" : "FAILED") << std::endl;

        // Test ray that starts inside the AABB
        Ray ray3(Vec3(0.5, 0.5, 0.5), Vec3(1, 0, 0)); // Ray starts inside the AABB
        float t3;
        bool intersect3 = aabb.rayAABB(ray3, t3);
        std::cout << "Test 3 (Ray starts inside AABB): " << (intersect3 ? "PASSED" : "FAILED") << std::endl;
    }

    // Test AABB surface area
    {
        AABB aabb;
        aabb.min = Vec3(0, 0, 0);
        aabb.max = Vec3(1, 1, 1);

        float area = aabb.area();
        std::cout << "Test 4 (AABB surface area): " << (std::abs(area - 6.0f) < 0.001f ? "PASSED" : "FAILED") << std::endl;
    }

    // Test Sphere ray intersection
    {
        Sphere sphere;
        sphere.init(Vec3(0, 0, 0), 1.0f); // Sphere at origin with radius 1

        // Test ray that intersects the sphere
        Ray ray1(Vec3(0, 0, -2), Vec3(0, 0, 1)); // Ray pointing towards the sphere
        float t1;
        bool intersect1 = sphere.rayIntersect(ray1, t1);
        std::cout << "Test 5 (Ray intersects sphere): " << (intersect1 ? "PASSED" : "FAILED") << std::endl;

        // Test ray that misses the sphere
        Ray ray2(Vec3(2, 2, 2), Vec3(0, 0, 1)); // Ray pointing away from the sphere
        float t2;
        bool intersect2 = sphere.rayIntersect(ray2, t2);
        std::cout << "Test 6 (Ray misses sphere): " << (!intersect2 ? "PASSED" : "FAILED") << std::endl;

        // Test ray that starts inside the sphere
        Ray ray3(Vec3(0, 0, 0), Vec3(0, 0, 1)); // Ray starts at the sphere's centre
        float t3;
        bool intersect3 = sphere.rayIntersect(ray3, t3);
        std::cout << "Test 7 (Ray starts inside sphere): " << (intersect3 ? "PASSED" : "FAILED") << std::endl;
    }

    // Test Plane ray intersection
    {
        Plane plane;
        Vec3 planeNormal(0, 1, 0); // Create a named Vec3 object
        plane.init(planeNormal, 0.0f); // Pass the named object

        // Test ray that intersects the plane
        Ray ray1(Vec3(0, 1, 0), Vec3(0, -1, 0)); // Ray pointing towards the plane
        float t1;
        bool intersect1 = plane.rayIntersect(ray1, t1);
        std::cout << "Test 8 (Ray intersects plane): " << (intersect1 ? "PASSED" : "FAILED") << std::endl;

        // Test ray that misses the plane (parallel)
        Ray ray2(Vec3(0, 1, 0), Vec3(1, 0, 0)); // Ray parallel to the plane
        float t2;
        bool intersect2 = plane.rayIntersect(ray2, t2);
        std::cout << "Test 9 (Ray misses plane): " << (!intersect2 ? "PASSED" : "FAILED") << std::endl;

        // Test ray that starts on the plane
        Ray ray3(Vec3(0, 0, 0), Vec3(0, 1, 0)); // Ray starts on the plane
        float t3;
        bool intersect3 = plane.rayIntersect(ray3, t3);
        std::cout << "Test 10 (Ray starts on plane): " << (intersect3 ? "PASSED" : "FAILED") << std::endl;
    }

    std::cout << "Tests completed." << std::endl;
}


int main(int argc, char* argv[])
{
    // Call tests
    runTests();

    // Initialize default parameters
    //std::string sceneName = "cornell-box";
    std::string sceneName = "Sibenik";
    std::string filename = "GI.hdr";
    unsigned int SPP = 8192;

    if (argc > 1)
    {
        std::unordered_map<std::string, std::string> args;
        for (int i = 1; i < argc; ++i)
        {
            std::string arg = argv[i];
            if (!arg.empty() && arg[0] == '-')
            {
                std::string argName = arg;
                if (i + 1 < argc)
                {
                    std::string argValue = argv[++i];
                    args[argName] = argValue;
                }
                else
                {
                    std::cerr << "Error: Missing value for argument '" << arg << "'\n";
                }
            }
            else
            {
                std::cerr << "Warning: Ignoring unexpected argument '" << arg << "'\n";
            }
        }
        for (const auto& pair : args)
        {
            if (pair.first == "-scene")
            {
                sceneName = pair.second;
            }
            if (pair.first == "-outputFilename")
            {
                filename = pair.second;
            }
            if (pair.first == "-SPP")
            {
                SPP = stoi(pair.second);
            }
        }
    }
    Scene* scene = loadScene(sceneName);
    GamesEngineeringBase::Window canvas;
    canvas.create((unsigned int)scene->camera.width, (unsigned int)scene->camera.height, "Tracer", 1.0f);
    RayTracer rt;
    rt.init(scene, &canvas);
    bool running = true;
    GamesEngineeringBase::Timer timer;
    while (running)
    {
        canvas.checkInput();
        canvas.clear();
        if (canvas.keyPressed(VK_ESCAPE))
        {
            break;
        }
        if (canvas.keyPressed('W'))
        {
            viewcamera.forward();
            rt.clear();
        }
        if (canvas.keyPressed('S'))
        {
            viewcamera.back();
            rt.clear();
        }
        if (canvas.keyPressed('A'))
        {
            viewcamera.left();
            rt.clear();
        }
        if (canvas.keyPressed('D'))
        {
            viewcamera.right();
            rt.clear();
        }
        if (canvas.keyPressed('E'))
        {
            viewcamera.flyUp();
            rt.clear();
        }
        if (canvas.keyPressed('Q'))
        {
            viewcamera.flyDown();
            rt.clear();
        }
        // Time how long a render call takes
        timer.reset();
        rt.render();
        float t = timer.dt();
        // Write
        std::cout << t << std::endl;
        if (canvas.keyPressed('P'))
        {
            rt.saveHDR(filename);
        }
        if (canvas.keyPressed('L'))
        {
            size_t pos = filename.find_last_of('.');
            std::string ldrFilename = filename.substr(0, pos) + ".png";
            rt.savePNG(ldrFilename);
        }
        if (SPP == rt.getSPP())
        {
            rt.saveHDR(filename);
            break;
        }
        canvas.present();
    }
    return 0;
}