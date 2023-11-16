#include "types.h"
#include <Eigen/Dense>
#include <vector>

#ifndef TIME_INTEGRATION
#define TIME_INTEGATION

namespace learnSPH
{
    namespace timeIntegration
    {
        class semiImplicitEuler
        {
            public:
              double radius;
              double v_max           = 1.0;
              bool boundary_checking = false;
              std::vector<types::boundary> boundaries;
              semiImplicitEuler(double radius);
              semiImplicitEuler(double radius, std::vector<types::boundary> boundaries);
              void integrationStep(std::vector<Eigen::Vector3d> &positions,
                                  std::vector<Eigen::Vector3d> &velocity,
                                  std::vector<Eigen::Vector3d> &accelerations, std::vector<bool> &deleteFlag,
                                  double dt, int &count_del);
          };
      } // namespace timeIntegration
} // namespace learnSPH

#endif  