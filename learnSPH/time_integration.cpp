#include "time_integration.h"

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

#include "types.h"

learnSPH::timeIntegration::semiImplicitEuler::semiImplicitEuler(double radius) {
  this->radius = radius;
  this->v_max = 0;
}

void learnSPH::timeIntegration::semiImplicitEuler::integrationStep(
    std::vector<Eigen::Vector3d> &positions,
    std::vector<Eigen::Vector3d> &velocity,
    std::vector<Eigen::Vector3d> &accelerations,
    double dt) 
{
    
  for (int i = 0; i < positions.size(); i++) {
    velocity[i] = velocity[i] + dt * (gravity + accelerations[i]);
    positions[i] = positions[i] + dt * velocity[i];

    if (velocity[i].norm() > v_max)
      v_max = velocity[i].norm();

    // if (positions[i].x() > x_boundary_max ||
    //     positions[i].x() < x_boundary_min ||
    //     positions[i].y() > y_boundary_max ||
    //     positions[i].y() < y_boundary_min ||
    //     positions[i].z() > z_boundary_max ||
    //     positions[i].z() < z_boundary_min) {
    //   positions.erase(positions.begin() + i);
    //   velocity.erase(velocity.begin() + i);
    //   accelerations.erase(accelerations.begin() + i);
    // }
  }
}