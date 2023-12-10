#include "theta_functions.h"
#include "../learnSPH/io.h"
#include "kernel.h"
#include <Eigen/Dense>
#include <cmath>
#include <cstdint>
#include <memory>
#include <sys/types.h>
#include <unordered_map>
#include <vector>

learnSPH::theta_functions::FluidThetaFunction::FluidThetaFunction(
    learnSPH::kernel::CubicSplineKernel &kernel, double c, double cell_width, double support_radius,
    uint nv_x, uint n_vy, uint n_vz)
    : kernel(kernel), c(c), cell_width(cell_width), support_radius(support_radius), n_vy(n_vy),
      n_vz(n_vz), n_vx(nv_x)
{
}

void learnSPH::theta_functions::FluidThetaFunction::computeLevelSet(
    std::vector<double> &level_set, std::vector<Eigen::Vector3d> &positions,
    std::vector<double> &densities, Eigen::Vector3d bound_min)
{
    uint vertexIdx;

    for (int pos = 0; pos < positions.size(); pos++) {

        int lower_x_abb =
            std::ceil((positions[pos].x() - support_radius - bound_min.x()) / cell_width);
        int upper_x_abb =
            std::floor((positions[pos].x() + support_radius - bound_min.x()) / cell_width);

        int lower_y_abb =
            std::ceil((positions[pos].y() - support_radius - bound_min.y()) / cell_width);
        int upper_y_abb =
            std::floor((positions[pos].y() + support_radius - bound_min.y()) / cell_width);

        int lower_z_abb =
            std::ceil((positions[pos].z() - support_radius - bound_min.z()) / cell_width);
        int upper_z_abb =
            std::floor((positions[pos].z() + support_radius - bound_min.z()) / cell_width);

        for (int i = lower_x_abb; i < upper_x_abb + 1; i++) {
            for (int j = lower_y_abb; j < upper_y_abb + 1; j++) {
                for (int k = lower_z_abb; k < upper_z_abb + 1; k++) {

                    Eigen::Vector3d vertex_pos =
                        Eigen::Vector3d(i * cell_width, j * cell_width, k * cell_width) + bound_min;

                    if ((positions[pos] - vertex_pos).squaredNorm() <
                        support_radius * support_radius) {
                        vertexIdx = i * n_vy * n_vz + j * n_vz + k;
                        level_set[vertexIdx] += (1 / densities[pos]) *
                                                kernel.kernel_function(positions[pos] - vertex_pos);
                    }
                }
            }
        }
    }
}

void learnSPH::theta_functions::FluidThetaFunction::computeLevelMap(
    std::unordered_map<uint64_t, double> &level_map, std::vector<Eigen::Vector3d> &positions,
    std::vector<double> &densities, Eigen::Vector3d bound_min)
{
    uint vertexIdx;

    for (int pos = 0; pos < positions.size(); pos++) {

        int lower_x_abb = std::floor(
            std::abs((positions[pos].x() - support_radius - bound_min.x()) / cell_width));
        int upper_x_abb =
            std::ceil(std::abs((positions[pos].x() + support_radius - bound_min.x()) / cell_width));

        int lower_y_abb = std::floor(
            std::abs((positions[pos].y() - support_radius - bound_min.y()) / cell_width));
        int upper_y_abb =
            std::ceil(std::abs((positions[pos].y() + support_radius - bound_min.y()) / cell_width));

        int lower_z_abb = std::floor(
            std::abs((positions[pos].z() - support_radius - bound_min.z()) / cell_width));
        int upper_z_abb =
            std::ceil(std::abs((positions[pos].z() + support_radius - bound_min.z()) / cell_width));

        for (int i = lower_x_abb; i < upper_x_abb + 1; i++) {
            for (int j = lower_y_abb; j < upper_y_abb + 1; j++) {
                for (int k = lower_z_abb; k < upper_z_abb + 1; k++) {

                    Eigen::Vector3d vertex_pos =
                        Eigen::Vector3d(i * cell_width, j * cell_width, k * cell_width) + bound_min;

                    if ((positions[pos] - vertex_pos).squaredNorm() <
                        support_radius * support_radius) {
                        vertexIdx = i * n_vy * n_vz + j * n_vz + k;
                        level_map[vertexIdx] += (1 / densities[pos]) *
                                                kernel.kernel_function(positions[pos] - vertex_pos);
                    }
                }
            }
        }
    }

    for (auto &idx : level_map) {
        level_map[idx.first] += -c;
    }
}