#pragma once
#include "../learnSPH/kernel.h"
#include <Eigen/Dense>
#include <cstdint>
#include <memory>
#include <sys/types.h>
#include <unordered_map>
#include <vector>

namespace learnSPH
{
namespace theta_functions
{
class FluidThetaFunction
{
  public:
    double cell_width;
    double c;
    double support_radius;
    std::vector<double> densities;
    uint n_vx, n_vy, n_vz;
    learnSPH::kernel::CubicSplineKernel &kernel;

    FluidThetaFunction(learnSPH::kernel::CubicSplineKernel &kernel, double c, double cell_width,
                       double support_radius, uint n_vx, uint n_vy, uint n_vz);
    void computeLevelSet(std::vector<double> &level_set, std::vector<Eigen::Vector3d> &positions,
                         std::vector<double> &densities, Eigen::Vector3d bound_min);
    void computeLevelMap(std::unordered_map<uint64_t, double> &level_map,
                         std::vector<Eigen::Vector3d> &positions, std::vector<double> &densities,
                         Eigen::Vector3d bound_min);
    double singleSignedDistance(Eigen::Vector3d pos, std::vector<Eigen::Vector3d> &positions,
                                std::vector<double> &densities, Eigen::Vector3d bound_min);
};
struct ssdArgs
{
    std::unique_ptr<std::vector<double>> densities;
};
} // namespace theta_functions
} // namespace learnSPH