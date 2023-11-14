#include "densities.h"
#include <iostream>

#include "../extern/CompactNSearch/include/CompactNSearch/CompactNSearch"
#include "../learnSPH/kernel.h"


void learnSPH::densities::compute_boundary_masses(
    std::vector<double>& output,
    std::vector<Eigen::Vector3d>& boundary_particles,
    unsigned int point_set_id,
    CompactNSearch::PointSet const& pointset,
    double density,
    learnSPH::kernel::CubicSplineKernel &cubic_kernel) 
{
    for (int i = 0; i < pointset.n_points(); ++i)
    {
        double kernel_sum = 0.0;
        kernel_sum += cubic_kernel.kernel_function(boundary_particles[i] - boundary_particles[i]);

        for (size_t j = 0; j < pointset.n_neighbors(point_set_id, i); ++j)
        {
            const unsigned int pid = pointset.neighbor(point_set_id, i, j);
            kernel_sum += cubic_kernel.kernel_function(boundary_particles[i] - boundary_particles[pid]);
        }
        output[i] = density / kernel_sum;
    }
}

void learnSPH::densities::compute_fluid_density(
    std::vector<double>& output,
    std::vector<Eigen::Vector3d>& particles,
    std::vector<Eigen::Vector3d>& boundary_particles,
    std::vector<double>& boundary_particles_masses,
    unsigned int point_set_id_fluid,
    CompactNSearch::PointSet const& ps_fluid,
    unsigned int point_set_id_boundary,
    CompactNSearch::PointSet const& ps_boundary,
    const double fluid_mass,
    learnSPH::kernel::CubicSplineKernel &cubic_kernel
    ) 
{
    for (int i = 0; i < ps_fluid.n_points(); ++i)
    {
        double density_sum = 0.0;

        density_sum += fluid_mass * cubic_kernel.kernel_function(particles[i] - particles[i]);

        // Get fluid neighbors of fluid point set.
        for (size_t j = 0; j < ps_fluid.n_neighbors(point_set_id_fluid, i); ++j)
        {
            const unsigned int pid = ps_fluid.neighbor(point_set_id_fluid, i, j);
            density_sum += fluid_mass * cubic_kernel.kernel_function(particles[i] - particles[pid]);
        }

        // Get boundary neighbors of fluid point set.
        for (size_t j = 0; j < ps_fluid.n_neighbors(point_set_id_boundary, i); ++j)
        {
            const unsigned int pid = ps_fluid.neighbor(point_set_id_boundary, i, j);
            density_sum += boundary_particles_masses[pid] * 
                            cubic_kernel.kernel_function(particles[i] - boundary_particles[pid]);
        }
        output[i] = density_sum;
    }
}