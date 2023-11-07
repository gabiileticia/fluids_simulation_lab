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
        double x_boundary_max, x_boundary_min, y_boundary_max, y_boundary_min, z_boundary_max, z_boundary_min;
        const Eigen::Vector3d gravity = {0, 0, -9.81};
        semiImplicitEuler(double radius);
        void integrationStep(std::vector<Eigen::Vector3d> &positions,
                            std::vector<Eigen::Vector3d> &velocity,
                            std::vector<Eigen::Vector3d> &accelerations,
                            double dt);
    };
  } // namespace timeIntegration
} // namespace learnSPH

#endif