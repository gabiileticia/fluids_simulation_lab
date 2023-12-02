#include "theta_functions.h"
#include "../learnSPH/io.h"
#include "kernel.h"
#include <Eigen/Dense>
#include <cmath>
#include <memory>
#include <sys/types.h>
#include <vector>

learnSPH::theta_functions::FluidThetaFunction::FluidThetaFunction(
    learnSPH::kernel::CubicSplineKernel &kernel, double c, double cell_width, double support_radius,
    uint n_vy, uint n_vz)
    : kernel(kernel), c(c), cell_width(cell_width), support_radius(support_radius), n_vy(n_vy),
      n_vz(n_vz)
{
}

void learnSPH::theta_functions::FluidThetaFunction::computeLevelMap(
    std::unordered_map<uint64_t, double> &level_map, std::vector<Eigen::Vector3d> &positions,
    std::vector<double> &densities, Eigen::Vector3d bound_min)
{
    uint vertexIdx;
    uint maxIdx = 0;
    uint hi, hj, hk;
    int track_lower, track_higher;
    Eigen::Vector3d track_pos;

    for (int pos = 0; pos < positions.size(); pos++) {
        int lower_x_abb =
            std::abs(std::ceil((positions[pos].x() - support_radius - bound_min.x()) / cell_width));
        int upper_x_abb =
            std::abs(std::floor((positions[pos].x() + support_radius - bound_min.x()) / cell_width));

        int lower_y_abb =
            std::abs(std::ceil((positions[pos].y() - support_radius - bound_min.y()) / cell_width));
        int upper_y_abb =
            std::abs(std::floor((positions[pos].y() + support_radius - bound_min.y()) / cell_width));

        int lower_z_abb =
            std::abs(std::ceil((positions[pos].z() - support_radius - bound_min.z()) / cell_width));
        int upper_z_abb =
            std::abs(std::floor((positions[pos].z() + support_radius - bound_min.z()) / cell_width));

        for (int i = lower_x_abb; i < upper_x_abb + 1; i++) {
            for (int j = lower_y_abb; j < upper_y_abb + 1; j++) {
                for (int k = lower_z_abb; k < upper_z_abb + 1; k++) {

                    Eigen::Vector3d vertex_pos = {i * cell_width + bound_min.x(),
                                                  j * cell_width + bound_min.y(),
                                                  k * cell_width + bound_min.z()};

                    if ((positions[pos] - vertex_pos).norm() < support_radius) {
                        vertexIdx = i * n_vy * n_vz + j * n_vz + k;
                        Eigen::Vector3d dx = positions[pos] - vertex_pos;
                        double kx          = kernel.kernel_function(dx);
                        level_map[vertexIdx] += (1 / densities[pos]) * kx;
                    }
                }
            }
        }
    }
    for(auto& idx : level_map){
        level_map[idx.first] += -c;
    }
    // std::cout << "Max found index: " << maxIdx << "\n"
    //           << "i: " << hi << "\nhj: " << hj << "\nhk" << hk << "\nbound_min: " << bound_min
    //           << "\nlower i: " << track_lower << "\nupper i: " << track_higher
    //           << "\nsupport radius: " << support_radius << "\ncell width: " << cell_width
    //           << "\npos: (" << track_pos.x() << ", " << track_pos.y() << ", " << track_pos.z()
    //           << ")\n";
}

void learnSPH::theta_functions::FluidThetaFunction::computeLevelSet(
    std::vector<double> &level_set, std::vector<Eigen::Vector3d> &positions,
    std::vector<double> &densities, Eigen::Vector3d bound_min)
{
    // dimensions
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

                    if ((positions[pos] - vertex_pos).norm() < support_radius) {
                        vertexIdx = i * n_vy * n_vz + j * n_vz + k;
                        level_set[vertexIdx] += (1 / densities[pos]) *
                                                kernel.kernel_function(positions[pos] - vertex_pos);
                    }
                }
            }
        }
    }
}

double learnSPH::theta_functions::FluidThetaFunction::singleSignedDistance(
    Eigen::Vector3d pos, std::vector<Eigen::Vector3d> &positions, std::vector<double> &densities,
    Eigen::Vector3d bound_min)
{
    uint vertexIdx;
    double signedDistance;

    int lower_x_abb = std::ceil(pos.x() - support_radius - bound_min.x() / cell_width);
    int upper_x_abb = std::ceil(pos.x() + support_radius - bound_min.x() / cell_width);
    int lower_y_abb = std::ceil(pos.y() - support_radius - bound_min.y() / cell_width);
    int upper_y_abb = std::ceil(pos.y() + support_radius - bound_min.y() / cell_width);
    int lower_z_abb = std::ceil(pos.z() - support_radius - bound_min.z() / cell_width);
    int upper_z_abb = std::ceil(pos.z() + support_radius - bound_min.z() / cell_width);

    for (int i = lower_x_abb; i < upper_x_abb; i++) {
        for (int j = lower_y_abb; j < upper_y_abb; j++) {
            for (int k = lower_z_abb; k < upper_z_abb; k++) {
                Eigen::Vector3d vertexPos = {i * cell_width + bound_min.x(),
                                             j * cell_width + bound_min.y(),
                                             k * cell_width + bound_min.z()};

                if ((pos - vertexPos).norm() < support_radius) {
                    vertexIdx = i * n_vy * n_vz + j * n_vz + k;
                    signedDistance +=
                        (1 / (densities)[i]) * kernel.kernel_function(pos - vertexPos);
                }
            }
        }
    }
    return signedDistance;
}

void learnSPH::theta_functions::Torus::computeLevelSet(std::vector<double> &level_set)
{
    for (int i = 0; i < this->nx; i++) {
        for (int j = 0; j < this->ny; j++) {
            for (int k = 0; k < this->nz; k++) {
                Eigen::Vector3d pos = {i * cellWidth + origin.x(), j * cellWidth + origin.y(),
                                       k * cellWidth + origin.z()};
                uint vertexIdx      = i * this->ny * this->nz + j * this->nz + k;
                double interval(std::sqrt(pos.x() * pos.x() + pos.y() * pos.y() - this->R));
                level_set[vertexIdx] = this->r * this->r - interval * interval - pos.z() * pos.z();
            }
        }
    }
}

void learnSPH::theta_functions::Torus::computeLevelMap(
    std::unordered_map<uint64_t, double> &level_map)
{
    for (int i = 0; i < this->nx; i++) {
        for (int j = 0; j < this->ny; j++) {
            for (int k = 0; k < this->nz; k++) {
                Eigen::Vector3d pos = {i * cellWidth + origin.x(), j * cellWidth + origin.y(),
                                       k * cellWidth + origin.z()};
                uint vertexIdx      = i * this->ny * this->nz + j * this->nz + k;
                double interval(std::sqrt(pos.x() * pos.x() + pos.y() * pos.y() - this->R));
                level_map[vertexIdx] = this->r * this->r - interval * interval - pos.z() * pos.z();
            }
        }
    }
}

double learnSPH::theta_functions::Torus::singleSignedDistance(Eigen::Vector3d pos)
{
    double interval = (std::sqrt(pos.x() * pos.x() + pos.y() * pos.y()) - this->R);
    return this->r * this->r - interval * interval - pos.z() * pos.z();
}