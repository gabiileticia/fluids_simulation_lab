#include "types.h"
#include <Eigen/Dense>
#include <array>
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
              std::vector<types::object> boundaries;
              Eigen::Vector3d simBoundary_min;
              Eigen::Vector3d simBoundary_max;

              semiImplicitEuler(double radius, std::vector<types::object> objects, std::vector<Eigen::Vector3d> simBoundary);
              void integrationStep(std::vector<Eigen::Vector3d> &positions,
                                  std::vector<Eigen::Vector3d> &velocity,
                                  std::vector<Eigen::Vector3d> &accelerations, std::vector<bool> &deleteFlag,
                                  double dt, int &count_del, Eigen::Vector3d &min_fluid_reco, Eigen::Vector3d &max_fluid_reco);
          };
      } // namespace timeIntegration
} // namespace learnSPH

#endif  