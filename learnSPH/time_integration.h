#include "types.h"
#include <Eigen/Dense>
#include <vector>

#ifndef TIME_INTEGRATION
#define TIME_INTEGATION

namespace learnSPH {
namespace timeIntegration {
class semiImplicitEuler {
public:
  double radius;
  double v_max = 1.0;
  bool boundary_checking = false;
  Eigen::Vector3d max_boundary;
  Eigen::Vector3d min_boundary;
  semiImplicitEuler(double radius);
  semiImplicitEuler(double radius, Eigen::Vector3d max_boundary,
                    Eigen::Vector3d min_boundary);
  void integrationStep(std::vector<Eigen::Vector3d> &positions,
                       std::vector<Eigen::Vector3d> &velocity,
                       std::vector<Eigen::Vector3d> &accelerations,
                       std::vector<double> &densities,
                       double dt);
};
} // namespace timeIntegration
} // namespace learnSPH

#endif