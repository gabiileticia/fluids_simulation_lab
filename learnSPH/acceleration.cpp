#include "acceleration.h"

#include "kernel.h"
#include <algorithm>
#include <vector>

learnSPH::acceleration::Acceleration::Acceleration(
    double B, double v_f, double v_b, double h, double roh_0,
    Eigen::Vector3d gravity, learnSPH::kernel::CubicSplineKernel &kernel)
    : B(B), v_f(v_f), v_b(v_b), h(h), kernel(kernel), roh_0(roh_0),
      gravity(gravity) {}

void learnSPH::acceleration::Acceleration::pressure(
    std::vector<double> &particles_pressure,
    std::vector<double> &particles_density, double rest_density) {
  double delta_density;
  for (int i = 0; i < particles_density.size(); i++) {
    delta_density = this->B * particles_density[i] - rest_density;
    particles_pressure[i] = (delta_density > 0.0) ? delta_density : 0.0;
  }
}

void learnSPH::acceleration::Acceleration::accelerations(
    std::vector<Eigen::Vector3d> &accelerations,
    std::vector<double> &roh, // density
    std::vector<double> &p,   // pressure
    unsigned int ps_id_fluid, unsigned int ps_id_boundary,
    CompactNSearch::PointSet &fluid_neighbors,
    CompactNSearch::PointSet &boundary_neighbors,
    std::vector<Eigen::Vector3d> &fluid_particles,
    std::vector<Eigen::Vector3d> &boundary_particles,
    std::vector<Eigen::Vector3d> &velocity,
    std::vector<double> &boundary_masses, const double fluid_mass) {
  Eigen::Vector3d ap, av, ae; // pressure, velocity and mass forces
  Eigen::Vector3d ff_inter_pressure, fs_inter_pressure, ff_inter_velocity,
      fs_inter_velocity; // fluid fluid and fluid solid interaction
  Eigen::Vector3d dx;
  ff_inter_pressure = {0, 0, 0};
  fs_inter_pressure = {0, 0, 0};
  ff_inter_velocity = {0, 0, 0};
  fs_inter_velocity = {0, 0, 0};

  for (int i = 0; i < fluid_particles.size(); ++i) {
    for (int j = 0; j < fluid_neighbors.n_neighbors(ps_id_fluid, i); j++) {
      unsigned int pid = fluid_neighbors.neighbor(ps_id_fluid, i, j);
      ff_inter_pressure +=
          fluid_mass *
          (p[i] / (roh[i] * roh[i]) + p[pid] / (roh[pid] * roh[pid])) *
          kernel.kernel_gradient(fluid_particles[i] - fluid_particles[pid]);

      dx = fluid_particles[i] - fluid_particles[pid];

      ff_inter_velocity += fluid_mass / roh[pid] *
                           (velocity[i] - velocity[pid]) * (dx).transpose() *
                           kernel.kernel_gradient(dx) /
                           (dx.norm() * dx.norm() + 0.01 * this->h * this->h);
    }
    for (int k = 0; k < boundary_neighbors.n_neighbors(ps_id_boundary, i);
         ++k) {
      unsigned int pid = boundary_neighbors.neighbor(ps_id_boundary, i, k);
      dx = fluid_particles[i] - boundary_particles[pid];
      fs_inter_pressure += this->roh_0 * boundary_masses[pid] *
                           (p[i] / roh[i] * roh[i]) *
                           kernel.kernel_gradient(dx);

      fs_inter_velocity += boundary_masses[pid] * velocity[i] * dx.transpose() *
                           kernel.kernel_gradient(dx) /
                           (dx.norm() * dx.norm() + 0.01 * h * h);
    }
    ap = -ff_inter_pressure - fs_inter_pressure;
    av = 2 * this->v_f * ff_inter_velocity + 2 * this->v_b * fs_inter_velocity;
    ae = this->gravity;

    accelerations[i] = ap + av + ae;
  }
}