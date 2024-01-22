#include "surface_tension.h"

void learnSPH::surface_tension::compute_forces(
                            std::vector<Eigen::Vector3d> &surface_tension_forces,
                            double gama, 
                            double adhesion_coefficient, 
                            const double fluid_mass,
                            const double rest_density,
                            std::vector<Eigen::Vector3d> &smoothed_color_field,
                            learnSPH::kernel::CubicSplineKernel &cubic_kernel,
                            std::vector<Eigen::Vector3d> &particles,
                            std::vector<Eigen::Vector3d> &boundary_particles,
                            std::vector<double> &particles_densities,
                            unsigned int point_set_id_fluid,
                            CompactNSearch::PointSet const &ps_fluid,
                            unsigned int point_set_id_boundary,
                            std::vector<double> boundary_mass)
{

    Eigen::Vector3d ff_cohesion, ff_curvature, surface_tension_sum, ff_adhesion; 
    double k_ij;

    surface_tension_forces.resize(particles_densities.size());
    for (int i = 0; i < ps_fluid.n_points(); ++i) {

        surface_tension_sum = {0,0,0};
        ff_adhesion = {0,0,0};
        int a =0;

        // Get fluid neighbors of fluid point set.
        for (size_t j = 0; j < ps_fluid.n_neighbors(point_set_id_fluid, i); ++j) {
            const unsigned int pid = ps_fluid.neighbor(point_set_id_fluid, i, j);

            ff_cohesion = - gama * fluid_mass 
                * cubic_kernel.cohesion_kernel_function((particles[i] - particles[pid]).norm())
                * (particles[i] - particles[pid]) 
                / ((particles[i] - particles[pid]).norm());
            
            ff_curvature = - gama * (smoothed_color_field[i] - smoothed_color_field[pid]);
            k_ij = 2 * rest_density / (particles_densities[i] + particles_densities[pid]);
            surface_tension_sum += k_ij * (ff_cohesion + ff_curvature);
        }

        // Get boundary neighbors of fluid point set.
        for (size_t j = 0; j < ps_fluid.n_neighbors(point_set_id_boundary, i); ++j) {
            const unsigned int pid = ps_fluid.neighbor(point_set_id_boundary, i, j);

            ff_adhesion -= adhesion_coefficient * fluid_mass * boundary_mass[pid] 
                * cubic_kernel.adhesion_kernel_function((particles[i] - boundary_particles[pid]).norm())
                * (particles[i] - boundary_particles[pid])
                / (particles[i] - boundary_particles[pid]).norm();
        }

        surface_tension_forces[i] = surface_tension_sum + ff_adhesion;
    }

}


void learnSPH::surface_tension::compute_smoothed_color_field(
            std::vector<Eigen::Vector3d> &smoothed_color_field,
            double c,
            const double fluid_mass,
            learnSPH::kernel::CubicSplineKernel &cubic_kernel,
            std::vector<double> &particles_densities,
            unsigned int point_set_id_fluid,
            CompactNSearch::PointSet const &ps_fluid,
            std::vector<Eigen::Vector3d> &particles
            )
{
    smoothed_color_field.resize(particles_densities.size());
    Eigen::Vector3d smoothed_color_sum; 
    for (int i = 0; i < ps_fluid.n_points(); ++i) {

        smoothed_color_sum = {0,0,0};

        // Get fluid neighbors of fluid point set.
        for (size_t j = 0; j < ps_fluid.n_neighbors(point_set_id_fluid, i); ++j) {
            const unsigned int pid = ps_fluid.neighbor(point_set_id_fluid, i, j);

            smoothed_color_sum += fluid_mass 
                * cubic_kernel.kernel_gradient(particles[i] - particles[pid]) 
                / particles_densities[pid];
        }
        smoothed_color_field[i] = c * smoothed_color_sum;
    }

}