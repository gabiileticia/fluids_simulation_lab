#include <Eigen/Dense>
#include <vector>
#include "types.h"

#ifndef TIME_INTEGRATION
#define TIME_INTEGATION

namespace learnSPH {
namespace timeIntegration {
class semiImplicitEuler {
public:
  double dt;
  Eigen::Vector3d vec_dt;
  double radius;
  double v_max;
  const double gravity = -9.81;
  semiImplicitEuler(double dt, double radius, double v_max);
  void integrationStep(std::vector<learnSPH::types::Particle> &particles);
};
} // namespace timeIntegration
} // namespace learnSPH

#endif