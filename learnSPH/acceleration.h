#include "../extern/CompactNSearch/include/CompactNSearch/CompactNSearch.h"
#include "kernel.h"
#include <Eigen/Dense>
#include <vector>

#ifndef ACCELERATION
#define ACCELERATION

namespace learnSPH {
namespace acceleration {
class Acceleration {
public:
  double B;     // stiffness constant
  double v_f;   // viscosity
  double v_b;   // friction with boundary
  double h;     // timestepsize
  double roh_0; // fluid rest density
  Eigen::Vector3d gravity;
  learnSPH::kernel::CubicSplineKernel kernel;

  Acceleration(double B, double v_f, double v_b, double h, double roh_0,
               Eigen::Vector3d gravity,
               learnSPH::kernel::CubicSplineKernel &kernel);

  void pressure(std::vector<double> &particles_pressure,
                std::vector<double> &particle_density, double rest_density);

  void accelerations(std::vector<Eigen::Vector3d> &accelerations,
                     std::vector<double> &particles_density,
                     std::vector<double> &particles_pressure,
                     unsigned int ps_id_fluid, unsigned int ps_is_boundary,
                     CompactNSearch::PointSet &fluid_neighbors,
                     CompactNSearch::PointSet &boundary_neighbors,
                     std::vector<Eigen::Vector3d> &fluid_particles,
                     std::vector<Eigen::Vector3d> &boundary_particles,
                     std::vector<Eigen::Vector3d> &velocity,
                     std::vector<double> &boundary_masses,
                     const double fluid_mass);
};
} // namespace acceleration
} // namespace learnSPH

#endif