#include "catch.hpp"
#include <Eigen/Dense>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "../learnSPH/io.h"
#include "../learnSPH/marching_cubes.h"
#include "../learnSPH/theta_functions.h"

TEST_CASE("Tests for marching cubes class", "[mcubes]")
{

    double cellWidth, epsilon, c;
    // dimensions
    uint nx, ny, nz;
    Eigen::Vector3d min, max;
    bool implicitFlag = false;

    double R  = 0.7;
    double r  = 0.2;
    cellWidth = 0.04;
    epsilon   = 1e-6;
    c         = -0.55;

    min = {-1, -1, -0.3};
    max = {1, 1, 0.3};
    // plus 1 so when it doesn't exactly fit we don't potentially loose parts of the shape
    nx = (max.x() - min.x()) / cellWidth + 1;
    ny = (max.y() - min.y()) / cellWidth + 1;
    nz = (max.z() - min.z()) / cellWidth + 1;

    std::vector<double> level_set(nx * ny * nz);
    std::unordered_map<uint64_t, double> levelMap(nx * ny * nz);

    learnSPH::theta_functions::Torus torus(min, r, R, cellWidth, nx, ny, nz);
    learnSPH::surface::MarchingCubes mcubes(cellWidth, nx, ny, nz, min, torus, epsilon,
                                            implicitFlag, c);

    torus.computeLevelSet(level_set);
    torus.computeLevelMap(levelMap);

    SECTION("Testing isosurface function")
    {
        mcubes.get_Isosurface(level_set);
        mcubes.get_Isosurface_sparse(levelMap);
    }
    SECTION("Testing Normals computation")
    {
        mcubes.get_Isosurface(level_set);
        mcubes.compute_normals();
    }

    SECTION("Testing output of file")
    {
        mcubes.get_Isosurface_sparse(levelMap);
        mcubes.compute_normals();
        learnSPH::write_tri_mesh_to_vtk("torus.vtk", mcubes.intersections, mcubes.triangles,
                                        mcubes.intersectionNormals);
        std::cout << "All finished";
    }
}