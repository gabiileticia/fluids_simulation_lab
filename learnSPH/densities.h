#pragma once
#include <array>
#include <vector>

#include <Eigen/Dense>
#include "../extern/CompactNSearch/include/CompactNSearch/CompactNSearch"
#include "kernel.h"


namespace learnSPH {
    namespace densities {
        void compute_boundary_masses(
            std::vector<double>& output,
            std::vector<Eigen::Vector3d>& boundary_particles,
            unsigned int point_set_id,
            CompactNSearch::PointSet const& pointset,
            const double density,
            learnSPH::kernel::CubicSplineKernel &kernel);
        void compute_fluid_density(
            std::vector<double>& output,
            std::vector<Eigen::Vector3d>& particles,
            std::vector<Eigen::Vector3d>& boundary_particles,
            std::vector<double>& boundary_particles_masses,
            unsigned int point_set_id_fluid,
            CompactNSearch::PointSet const& ps_fluid,
            unsigned int point_set_id_boundary,
            CompactNSearch::PointSet const& ps_boundary,
            const double fluid_mass,
            learnSPH::kernel::CubicSplineKernel &kernel
            );
    }
}