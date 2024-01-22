#include "acceleration.h"
#include "kernel.h"
#include <algorithm>
#include <cstdio>
#include <omp.h>
#include <stdlib.h> // rand
#include <vector>

learnSPH::acceleration::Acceleration::Acceleration(
    double B, double v_f, double v_b, double h, double roh_0, Eigen::Vector3d gravity,
    learnSPH::kernel::CubicSplineKernel &kernel, unsigned int ps_id_fluid,
    unsigned int ps_id_boundary, CompactNSearch::PointSet const& fluid_neighbors,
    CompactNSearch::PointSet const& boundary_neighbors, std::vector<Eigen::Vector3d> &fluid_particles,
    std::vector<Eigen::Vector3d> &accelerations, std::vector<Eigen::Vector3d> &boundary_particles,
    std::vector<Eigen::Vector3d> &velocity, std::vector<double> &boundary_mass,
    std::vector<double> &particle_density, std::vector<double> &particles_pressure)
    : B(B), v_f(v_f), v_b(v_b), h(h), kernel(kernel), roh_0(roh_0), gravity(gravity),
      ps_id_fluid(ps_id_fluid), ps_id_boundary(ps_id_boundary), fluid_neighbors(fluid_neighbors),
      boundary_neighbors(boundary_neighbors), accelerations(accelerations),
      boundary_particles(boundary_particles), fluid_particles(fluid_particles), velocity(velocity),
      roh(particle_density), p(particles_pressure), boundary_mass(boundary_mass)
{
    this->h_square = h * h;
}

void learnSPH::acceleration::Acceleration::pressure(double rest_density)
{
    double delta_density;

    for (int i = 0; i < roh.size(); i++) {
        delta_density = this->B * (roh[i] - rest_density);
        p[i]          = std::max(0.0, delta_density);
    }
}

void learnSPH::acceleration::Acceleration::wcsph_accelerations(double boundary_density,
                                                               const double fluid_particle_mass)
{

    Eigen::Vector3d ap, av, ae; // pressure, velocity and mass forces
    Eigen::Vector3d ff_inter_pressure, fs_inter_pressure, ff_inter_velocity,
        fs_inter_velocity; // fluid fluid and fluid solid interaction
    Eigen::Vector3d dx;
    Eigen::Vector3d kernel_grad;
    double dx_norm;
    double h_square_cent = 0.01 * this->h_square;
    double roh_square_i_inverse;
    double vol_boundary;

    for (int i = 0; i < fluid_neighbors.n_points(); ++i) {

        ff_inter_pressure    = {0, 0, 0};
        fs_inter_pressure    = {0, 0, 0};
        ff_inter_velocity    = {0, 0, 0};
        fs_inter_velocity    = {0, 0, 0};
        roh_square_i_inverse = 1 / (roh[i] * roh[i]);

        for (size_t j = 0; j < fluid_neighbors.n_neighbors(ps_id_fluid, i); j++) {

            unsigned int pid = fluid_neighbors.neighbor(ps_id_fluid, i, j);
            dx               = fluid_particles[i] - fluid_particles[pid];
            kernel_grad      = kernel.kernel_gradient(dx);
            dx_norm          = dx.norm();

            ff_inter_pressure += fluid_particle_mass *
                                 (p[i] * roh_square_i_inverse + p[pid] / (roh[pid] * roh[pid])) *
                                 kernel_grad;

            ff_inter_velocity += fluid_particle_mass * (velocity[i] - velocity[pid]) *
                                 (dx).transpose() * kernel_grad /
                                 (roh[pid] * (dx_norm * dx_norm + h_square_cent));
        }

        for (size_t j = 0; j < fluid_neighbors.n_neighbors(ps_id_boundary, i); ++j) {
            unsigned int pid = fluid_neighbors.neighbor(ps_id_boundary, i, j);
            dx               = fluid_particles[i] - boundary_particles[pid];
            kernel_grad      = kernel.kernel_gradient(dx);
            dx_norm          = dx.norm();
            vol_boundary     = boundary_mass[pid] / boundary_density;

            fs_inter_pressure +=
                this->roh_0 * vol_boundary * p[i] * roh_square_i_inverse * kernel_grad;

            fs_inter_velocity += vol_boundary * velocity[i] * dx.transpose() * kernel_grad /
                                 (dx_norm * dx_norm + h_square_cent);
        }

        ap = (-ff_inter_pressure - fs_inter_pressure);
        av = 2 * this->v_f * ff_inter_velocity + 2 * this->v_b * fs_inter_velocity;
        ae = this->gravity;

        accelerations[i] = ap + av + ae;
    }
}

void learnSPH::acceleration::Acceleration::pbf_accelerations(double boundary_density,
                                                             const double fluid_particle_mass)
{

    Eigen::Vector3d av, ae;                               // pressure, velocity and mass forces
    Eigen::Vector3d ff_inter_velocity, fs_inter_velocity; // fluid fluid and fluid solid interaction
    Eigen::Vector3d dx;
    Eigen::Vector3d kernel_grad;
    double dx_norm;
    double h_square_cent = 0.01 * this->h_square;
    double vol_boundary;

    for (int i = 0; i < fluid_neighbors.n_points(); ++i) {

        ff_inter_velocity = {0, 0, 0};
        fs_inter_velocity = {0, 0, 0};

        for (size_t j = 0; j < fluid_neighbors.n_neighbors(ps_id_fluid, i); j++) {

            unsigned int pid = fluid_neighbors.neighbor(ps_id_fluid, i, j);
            dx               = fluid_particles[i] - fluid_particles[pid];
            kernel_grad      = kernel.kernel_gradient(dx);
            dx_norm          = dx.norm();

            ff_inter_velocity += fluid_particle_mass * (velocity[i] - velocity[pid]) *
                                 (dx).transpose() * kernel_grad /
                                 (roh[pid] * (dx_norm * dx_norm + h_square_cent));
        }

        for (size_t j = 0; j < fluid_neighbors.n_neighbors(ps_id_boundary, i); ++j) {
            unsigned int pid = fluid_neighbors.neighbor(ps_id_boundary, i, j);
            dx               = fluid_particles[i] - boundary_particles[pid];
            kernel_grad      = kernel.kernel_gradient(dx);
            dx_norm          = dx.norm();
            vol_boundary     = boundary_mass[pid] / boundary_density;

            fs_inter_velocity += vol_boundary * velocity[i] * dx.transpose() * kernel_grad /
                                 (dx_norm * dx_norm + h_square_cent);
        }

        av = 2 * this->v_f * ff_inter_velocity + 2 * this->v_b * fs_inter_velocity;
        ae = this->gravity;

        accelerations[i] = av + ae;
    }
}