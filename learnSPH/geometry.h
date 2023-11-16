#pragma once
#include <array>
#include <vector>

#include <Eigen/Dense>

namespace learnSPH {
    namespace geometry {
        void load_n_sample_boundary(std::vector<Eigen::Vector3d>& output
                            , std::vector<std::string> file_name
                            , double boundary_sampling_distance);
        void load_n_sample_fluids(std::vector<Eigen::Vector3d>& output_positions
                            , std::vector<Eigen::Vector3d>& output_velocities
                            , std::vector<Eigen::Vector3d> fluid_begin
                            , std::vector<Eigen::Vector3d> fluid_end
                            , double fluid_sampling_distance
                            , std::vector<Eigen::Vector3d> fluid_velocities);
    }
}