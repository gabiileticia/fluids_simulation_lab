#pragma once
#include <array>
#include <vector>

#include <Eigen/Dense>

namespace learnSPH {
    namespace geometry {
        void load_n_sample_boundary(std::vector<Eigen::Vector3d>& output
                            , std::string file_name
                            , double boundary_sampling_distance);
    }
}