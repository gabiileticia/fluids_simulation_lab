#include "theta_functions.h"
#include <Eigen/Dense>
#include <vector>
#include "../learnSPH/io.h"
#include "kernel.h"

struct Radii
{
    double r;
    double R;
};

double learnSPH::theta_functions::ImplicitTorus(Eigen::Vector3d pos, void *args)
{
    // pointer conversion
    Radii *radii = ((Radii *)args);

    double interval = (std::sqrt(pos.x() * pos.x() + pos.y() * pos.y()) - radii->R);
    return radii->r * radii->r - interval * interval - pos.z() * pos.z();
}


std::vector<double> learnSPH::theta_functions::fluid_reconstruction(std::vector<Eigen::Vector3d> positions
, std::vector<double>  densities, Eigen::Vector3d bound_min, Eigen::Vector3d bound_max)
{
    double cellWidth, epsilon;
    // dimensions
    uint nx, ny, nz,n_vx, n_vy, n_vz,vertexIdx;

    learnSPH::kernel::CubicSplineKernel cubic_kernel(1.2*2.0*0.005);
    
    cellWidth = 0.005 * 1;
    // plus 1 so when it doesn't exactly fit we don't potentially loose parts of the shape
    nx = (bound_max.x() - bound_min.x()) / cellWidth + 1;
    ny = (bound_max.y() - bound_min.y()) / cellWidth + 1;
    nz = (bound_max.z() - bound_min.z()) / cellWidth + 1;

    n_vx = nx + 1;
    n_vy = ny + 1;
    n_vz = nz + 1;

    double c = 0.55;
    double support_radius = 0.024;
    std::vector<double> levelSet(n_vx * n_vy * n_vz, -c);

    for(int pos=0; pos < positions.size(); pos++){

        int lower_x_abb = std::ceil((positions[pos].x() - support_radius - bound_min.x())/cellWidth);
        int upper_x_abb = std::floor((positions[pos].x() + support_radius - bound_min.x())/cellWidth);

        int lower_y_abb = std::ceil((positions[pos].y() - support_radius - bound_min.y())/cellWidth);
        int upper_y_abb = std::floor((positions[pos].y() + support_radius - bound_min.y())/cellWidth);

        int lower_z_abb = std::ceil((positions[pos].z() - support_radius - bound_min.z())/cellWidth);
        int upper_z_abb = std::floor((positions[pos].z() + support_radius - bound_min.z())/cellWidth);

        for (int i=lower_x_abb; i <upper_x_abb + 1; i++){
            for (int j=lower_y_abb; j < upper_y_abb + 1; j++){
                for(int k=lower_z_abb; k < upper_z_abb + 1; k++){
                    
                    Eigen::Vector3d vertex_pos = Eigen::Vector3d(i*cellWidth,j*cellWidth,k*cellWidth);

                    if ((positions[pos] - vertex_pos).norm() < support_radius){
                        vertexIdx = i * n_vy * n_vz + j * n_vz + k;
                        levelSet[vertexIdx] += (1/densities[i]) * cubic_kernel.kernel_function(positions[pos] - vertex_pos);
                    }
                }
            }
        }
    }

    return levelSet;
}

