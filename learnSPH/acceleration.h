#include "../extern/CompactNSearch/include/CompactNSearch/CompactNSearch.h"
#include "CompactNSearch/PointSet.h"
#include "kernel.h"
#include <Eigen/Dense>
#include <vector>

#ifndef ACCELERATION
#define ACCELERATION

namespace learnSPH
{
namespace acceleration
{
class Acceleration
{
  public:
    double B;     // stiffness constant
    double v_f;   // viscosity
    double v_b;   // friction with boundary
    double h;     // timestepsize
    double roh_0; // fluid rest density
    double h_square;
    Eigen::Vector3d gravity;
    learnSPH::kernel::CubicSplineKernel kernel;
    unsigned int ps_id_fluid, ps_id_boundary;
    CompactNSearch::PointSet const &fluid_neighbors, &boundary_neighbors;
    std::vector<Eigen::Vector3d> &accelerations, &fluid_particles, &boundary_particles, &velocity;
    std::vector<double> &roh, &p, &boundary_mass;

    Acceleration(double B, double v_f, double v_b, double h, double roh_0, Eigen::Vector3d gravity,
                 learnSPH::kernel::CubicSplineKernel &kernel, unsigned int ps_id_fluid,
                 unsigned int ps_id_boundary, CompactNSearch::PointSet const& fluid_neighbors,
                 CompactNSearch::PointSet const& boundary_neighbors,
                 std::vector<Eigen::Vector3d> &fluid_particles,
                 std::vector<Eigen::Vector3d> &accelerations,
                 std::vector<Eigen::Vector3d> &boundary_particles,
                 std::vector<Eigen::Vector3d> &velocity,std::vector<double> &boundary_mass,
                 std::vector<double> &particle_density, std::vector<double> &particles_pressure);

    void pressure(double rest_density);

    void wcsph_accelerations(double boundary_density, const double fluid_mass);

    void pbf_accelerations(double boundary_density, const double fluid_particle_mass);
};
} // namespace acceleration
} // namespace learnSPH

#endif