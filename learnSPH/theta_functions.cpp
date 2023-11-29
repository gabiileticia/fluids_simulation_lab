#include "theta_functions.h"
#include <Eigen/Dense>
#include <vector>
#include "../learnSPH/io.h"
#include "kernel.h"


learnSPH::theta_functions::FluidThetaFunction::FluidThetaFunction(double c, double cell_width, double support_radius, uint n_vx, uint n_vy, uint n_vz){
    this->c = c;
    this->cell_width = cell_width;
    this->support_radius = support_radius;
    this->n_vx = n_vx;
    this->n_vy = n_vy;
    this->n_vz = n_vz;
}

void learnSPH::theta_functions::FluidThetaFunction::compute_fluid_reconstruction(std::vector<double> &level_set
, std::vector<Eigen::Vector3d> positions, std::vector<double>  densities, Eigen::Vector3d bound_min
, Eigen::Vector3d bound_max, learnSPH::kernel::CubicSplineKernel &kernel)
{
    // dimensions
    uint vertexIdx;

    for(int pos=0; pos < positions.size(); pos++){

        uint zero = 0;

        uint lower_x_abb = std::ceil((positions[pos].x() - support_radius - bound_min.x())/cell_width);
        uint upper_x_abb = std::floor((positions[pos].x() + support_radius - bound_min.x())/cell_width);

        uint lower_y_abb = std::ceil((positions[pos].y() - support_radius - bound_min.y())/cell_width);
        uint upper_y_abb = std::floor((positions[pos].y() + support_radius - bound_min.y())/cell_width);

        uint lower_z_abb = std::ceil((positions[pos].z() - support_radius - bound_min.z())/cell_width);
        uint upper_z_abb = std::floor((positions[pos].z() + support_radius - bound_min.z())/cell_width);

        for (int i=lower_x_abb; i < upper_x_abb + 1; i++){
            for (int j=lower_y_abb; j < upper_y_abb + 1; j++){
                for(int k=lower_z_abb; k < upper_z_abb + 1; k++){
                    
                    Eigen::Vector3d vertex_pos = Eigen::Vector3d(i*cell_width,j*cell_width,k*cell_width) + bound_min;

                    if ((positions[pos] - vertex_pos).norm() < support_radius){
                        vertexIdx = i * n_vy * n_vz + j * n_vz + k;
                        level_set[vertexIdx] += (1/densities[i]) * kernel.kernel_function(positions[pos] - vertex_pos);
                    }
                }
            }
        }
    }
}

