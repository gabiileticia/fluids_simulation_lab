#ifndef UTILS
#define UTILS

#include <Eigen/Dense>
#include <chrono>
#include <vector>

namespace learnSPH
{
namespace utils
{
void deleteOutOfBounds(std::vector<Eigen::Vector3d> &positions,
                       std::vector<Eigen::Vector3d> &velocity,
                       std::vector<Eigen::Vector3d> &accelerations, std::vector<double> &densities,
                       std::vector<double> &pressure, std::vector<bool> &deleteFlat,
                       int &count_del);
void create_simulation_folder(const std::string assign_number, std::string &timestamp);
void updateProgressBar(int &currentStep, int &maxSteps);;
} // namespace utils
} // namespace learnSPH

#endif