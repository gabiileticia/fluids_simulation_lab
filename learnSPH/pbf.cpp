#include "pbf.h"
#include "../learnSPH/kernel.h"
#include "../extern/CompactNSearch/include/CompactNSearch/CompactNSearch"

#include <vector>
#include <Eigen/Dense>

void learnSPH::pbf::compute_c(std::vector<double> &c, 
            std::vector<double> &particles_densities, 
            double rest_density)
{
    c.resize(particles_densities.size());
    for(int i = 0; i < c.size(); i++){
        c[i] = (particles_densities[i] / rest_density) - 1;
    }
}


void learnSPH::pbf::compute_s(std::vector<double> &s,
    std::vector<Eigen::Vector3d> &particles,
    std::vector<Eigen::Vector3d> &boundary_particles,
    const double fluid_mass,
    double rest_density,
    std::vector<double> boundary_mass,
    double boundary_density,
    learnSPH::kernel::CubicSplineKernel &cubic_kernel,
    unsigned int point_set_id_fluid,
    CompactNSearch::PointSet const &ps_fluid,
    unsigned int point_set_id_boundary)
{
    s.resize(particles.size());

    Eigen::Vector3d k_equals_i;
    double k_not_i;

    for (int i = 0; i < ps_fluid.n_points(); ++i) {
        
        k_equals_i = {0,0,0};
        k_not_i    = 0.0;

        // Get fluid neighbors of fluid point set.
        for (size_t j = 0; j < ps_fluid.n_neighbors(point_set_id_fluid, i); ++j) {
            const unsigned int pid = ps_fluid.neighbor(point_set_id_fluid, i, j);

            k_equals_i += (fluid_mass / rest_density) 
                * cubic_kernel.kernel_gradient(particles[i] - particles[pid]);

            k_not_i += ((- fluid_mass / rest_density) * cubic_kernel.kernel_gradient(particles[i] - particles[pid])).squaredNorm();
        }

        // Get boundary neighbors of fluid point set.
        for (size_t j = 0; j < ps_fluid.n_neighbors(point_set_id_boundary, i); ++j) {
            const unsigned int pid = ps_fluid.neighbor(point_set_id_boundary, i, j);

            k_equals_i += (boundary_mass[pid] / boundary_density) 
                * cubic_kernel.kernel_gradient(particles[i] - boundary_particles[pid]);
        }

        s[i] = (k_equals_i.squaredNorm() / fluid_mass) + (k_not_i / fluid_mass);
    }
}


void learnSPH::pbf::compute_lambda(
            std::vector<double> &lambda,
            std::vector<double> &c,
            std::vector<double> &s,
            double epsilon)
{
    lambda.resize(c.size());
    for (int i = 0 ; i < c.size(); i++){
        if (c[i] > 0)
            lambda[i] = - c[i] / (s[i] + epsilon);
        else
            lambda[i] = 0.0;
    }
}


void learnSPH::pbf::compute_dx(
            std::vector<Eigen::Vector3d> &dx,
            double rest_density,
            const double fluid_mass,
            std::vector<double> &lambda,
            learnSPH::kernel::CubicSplineKernel &cubic_kernel,
            std::vector<double> boundary_mass,
            unsigned int point_set_id_fluid,
            CompactNSearch::PointSet const &ps_fluid,
            unsigned int point_set_id_boundary,
            std::vector<Eigen::Vector3d> &particles,
            std::vector<Eigen::Vector3d> &boundary_particles)
{
    dx.resize(lambda.size());

    Eigen::Vector3d fluid_contribution;
    Eigen::Vector3d boundary_contribution;

    for (int i = 0; i < ps_fluid.n_points(); ++i) {

        fluid_contribution = {0,0,0};
        boundary_contribution = {0,0,0};
        
        // Get fluid neighbors of fluid point set.
        for (size_t j = 0; j < ps_fluid.n_neighbors(point_set_id_fluid, i); ++j) {
            const unsigned int pid = ps_fluid.neighbor(point_set_id_fluid, i, j);
            fluid_contribution += (lambda[i] + lambda[pid]) 
                * cubic_kernel.kernel_gradient(particles[i] - particles[pid]);
        }

        // Get boundary neighbors of fluid point set.
        for (size_t j = 0; j < ps_fluid.n_neighbors(point_set_id_boundary, i); ++j) {
            const unsigned int pid = ps_fluid.neighbor(point_set_id_boundary, i, j);
            boundary_contribution += (boundary_mass[pid] / fluid_mass) * lambda[i] 
                * cubic_kernel.kernel_gradient(particles[i] - boundary_particles[pid]);
        }

        dx[i] = (fluid_contribution + boundary_contribution) / rest_density;

    }

}


void learnSPH::pbf::update_pbf_positions(
            std::vector<Eigen::Vector3d> &particles,
            std::vector<Eigen::Vector3d> &dx)
{
    for(int i = 0; i < particles.size(); i++){
        particles[i] += dx[i];
    }
}


void learnSPH::pbf::update_positions_and_velocities(
            std::vector<Eigen::Vector3d> &particles,
            std::vector<Eigen::Vector3d> &pbf_particles_positions,
            std::vector<Eigen::Vector3d> &particles_velocities,
            double dt)
{
    for (int i = 0; i < particles.size(); i++){
        particles_velocities[i] = (pbf_particles_positions[i] - particles[i]) / dt;
        particles[i] = pbf_particles_positions[i];
    }
}