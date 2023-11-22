#include "catch.hpp"
#include <Eigen/Dense>

#include "../learnSPH/io.h"
#include "../learnSPH/marching_cubes.h"

struct Radii
{
    double r;
    double R;
};

double ImplicitTorus(Eigen::Vector3d pos, void *args)
{
    // pointer conversion
    Radii *radii = ((Radii *)args);

    double interval = (std::sqrt(pos.x() * pos.x() + pos.y() * pos.y()) - radii->R);
    return radii->r * radii->r - interval * interval - pos.z() * pos.z();
}

TEST_CASE("Tests for marching cubes class", "[mcubes]")
{

    double cellWidth, epsilon;
    // dimensions
    uint nx, ny, nz;
    Eigen::Vector3d min, max;
    Radii radii;

    radii.R   = 0.7;
    radii.r   = 0.2;
    cellWidth = 0.04;
    epsilon   = 1e-6;

    min = {-1, -1, -0.3};
    max = {1, 1, 0.3};
    // plus 1 so when it doesn't exactly fit we don't potentially loose parts of the shape
    nx = (max.x() - min.x()) / cellWidth + 1;
    ny = (max.y() - min.y()) / cellWidth + 1;
    nz = (max.z() - min.z()) / cellWidth + 1;

    learnSPH::surface::MarchingCubes mcubes(cellWidth, nx, ny, nz, min, ImplicitTorus, &radii,
                                            epsilon);

    SECTION("Testing isosurface function"){
        mcubes.get_Isosurface();
    }
    SECTION("Testing Normals computation"){
        mcubes.compute_normals();
    }

    SECTION("Testing output of file"){
        learnSPH::write_tri_mesh_to_vtk("torus.vtk", mcubes.intersections, mcubes.triangles, mcubes.intersectionNormals);
        std::cout << "All finished";
    }

}