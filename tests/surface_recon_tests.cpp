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
    uint vertexIdx;
    Eigen::Vector3d min, max;
    bool implicitFlag = true;
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

    std::vector<double> level_set((nx + 1) * (ny + 1) * (nz + 1));

    for(int i=0; i<nx + 1; i++){
        for(int j=0; j<ny + 1; j++){
            for (int k=0; k<nz + 1;k++){
                Eigen::Vector3d vertex_pos = Eigen::Vector3d(i*cellWidth,j*cellWidth,k*cellWidth) + min;

                vertexIdx = i * (ny + 1) * (nz + 1) + j * (nz + 1) + k;
                double circle = std::sqrt(vertex_pos.x() * vertex_pos.x() + vertex_pos.y() * vertex_pos.y()) - radii.R;
                level_set[vertexIdx] = radii.r * radii.r - circle * circle - vertex_pos.z() * vertex_pos.z();
                // if(level_set[vertexIdx] > 0){
                //     std::cout << vertex_pos.transpose() << std::endl;
                //     std::cout << level_set[vertexIdx]<< std::endl;
                // }
            }
        }
    }

    learnSPH::surface::MarchingCubes mcubes(cellWidth, nx, ny, nz, min,
                                            epsilon, implicitFlag);

    // SECTION("Testing isosurface function"){
    //     mcubes.get_isosurface(level_set);
    // }
    // SECTION("Testing Normals computation"){
    //     mcubes.get_isosurface(level_set);
    //     mcubes.compute_normals_gl(level_set);
    // }

    SECTION("Testing output of file"){
        mcubes.get_isosurface(level_set);
        std::cout << "iso" << std::endl;
        mcubes.compute_normals_gl(level_set);
        std::cout << "normals" << std::endl;
        learnSPH::write_tri_mesh_to_vtk("torus.vtk", mcubes.intersections, mcubes.triangles, mcubes.intersectionNormals);
        std::cout << "All finished";
    }

}