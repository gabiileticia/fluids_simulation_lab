#include <vector>
#include <Eigen/Dense>

#include "../learnSPH/kernel.h"
#include "../extern/CompactNSearch/include/CompactNSearch/CompactNSearch"

namespace learnSPH {
    namespace pbf {
        void compute_s(std::vector<double> &s,
            std::vector<Eigen::Vector3d> &particles,
            std::vector<Eigen::Vector3d> &boundary_particles,
            const double fluid_mass,
            double rest_density,
            std::vector<double> boundary_mass,
            double boundary_density,
            learnSPH::kernel::CubicSplineKernel &cubic_kernel,
            unsigned int point_set_id_fluid,
            CompactNSearch::PointSet const &ps_fluid,
            unsigned int point_set_id_boundary,
            CompactNSearch::PointSet const &ps_boundary
        );

        void compute_c(std::vector<double> &c, 
            std::vector<double> &particles_densities, 
            double rest_density
        );

        void compute_lambda(
            std::vector<double> &lambda,
            std::vector<double> &c,
            std::vector<double> &s,
            double epsilon
        );

        void compute_dx(
            std::vector<Eigen::Vector3d> &dx,
            double rest_density,
            const double fluid_mass,
            std::vector<double> &lambda,
            learnSPH::kernel::CubicSplineKernel &cubic_kernel,
            std::vector<double> boundary_mass,
            unsigned int point_set_id_fluid,
            CompactNSearch::PointSet const &ps_fluid,
            unsigned int point_set_id_boundary,
            CompactNSearch::PointSet const &ps_boundary,
            std::vector<Eigen::Vector3d> &particles,
            std::vector<Eigen::Vector3d> &boundary_particles
        );

        void update_pbf_positions(
            std::vector<Eigen::Vector3d> &particles,
            std::vector<Eigen::Vector3d> &dx
        );

        void update_positions_and_velocities(
            std::vector<Eigen::Vector3d> &particles,
            std::vector<Eigen::Vector3d> &pbf_particles_positions,
            std::vector<Eigen::Vector3d> &particles_velocities,
            double dt
        );
    }
}