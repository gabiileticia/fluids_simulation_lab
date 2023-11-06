#include "types.h"
#include <Eigen/Dense>
#include <vector>

#ifndef TIME_INTEGRATION
#define TIME_INTEGATION

namespace learnSPH {
namespace timeIntegration {
class semiImplicitEuler {
public:
  double dt;
  double radius;
  double v_max;
  const Eigen::Vector3d gravity = {0, 0, -9.81};
  semiImplicitEuler(double dt, double radius);
  void integrationStep(std::vector<Eigen::Vector3d> &positions,
                       std::vector<Eigen::Vector3d> &velocity,
                       std::vector<Eigen::Vector3d> &accelerations);
};
} // namespace timeIntegration
} // namespace learnSPH

#endif