#include "../learnSPH/kernel.h"
#include "../extern/CompactNSearch/include/CompactNSearch/CompactNSearch"

namespace learnSPH {
    namespace surface_tension {
        void compute_forces(
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
            std::vector<double> boundary_mass
        );
        void compute_smoothed_color_field(
            std::vector<Eigen::Vector3d> &smoothed_color_field,
            double c,
            const double fluid_mass,
            learnSPH::kernel::CubicSplineKernel &cubic_kernel,
            std::vector<double> &particles_densities,
            unsigned int point_set_id_fluid,
            CompactNSearch::PointSet const &ps_fluid,
            std::vector<Eigen::Vector3d> &particles
        );
    }
}