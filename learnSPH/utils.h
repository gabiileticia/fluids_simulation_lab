#ifndef UTILS
#define UTILS

#include <vector>
#include <Eigen/Dense>

namespace learnSPH {
namespace utils {
void deleteOutOfBounds(std::vector<Eigen::Vector3d> &positions,
std::vector<Eigen::Vector3d> &velocity,
std::vector<Eigen::Vector3d> &accelerations,
std::vector<double> &densities,
std::vector<double> &pressure,
std::vector<bool> &deleteFlat);
}
}



#endif