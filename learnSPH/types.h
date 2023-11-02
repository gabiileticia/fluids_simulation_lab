#include <Eigen/Dense>
#ifndef TYPES
#define TYPES

namespace learnSPH {
namespace types {
struct Particle {
  Eigen::Vector3d pos;
  Eigen::Vector3d velocity;
  Eigen::Vector3d force;
};
}
}

#endif