#pragma once
#include <array>
#include <vector>

#include <Eigen/Dense>

namespace learnSPH {
    namespace neighborhood_search {
        void brute_force_neighbors(const double h, std::vector<Eigen::Vector3d> particles, double beta);
        std::vector<std::string> brute_force_neighbor(const double h, std::vector<Eigen::Vector3d> particles, int particle_index, double beta);
        std::vector<std::vector<int>> brute_force_search(const double h, std::vector<Eigen::Vector3d> particles, double beta);
    }
}