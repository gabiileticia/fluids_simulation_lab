#include "time_integration.h"

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

#include "types.h"

learnSPH::timeIntegration::semiImplicitEuler::semiImplicitEuler(double dt,
                                                                double radius) {
  assert(dt <=
         0.5 * radius / std::pow(std::pow(v_max, 2) + std::pow(v_max, 2), 0.5));
  this->dt = dt;
  this->radius = radius;
  this->v_max = 0;
}

void learnSPH::timeIntegration::semiImplicitEuler::integrationStep(
    std::vector<Eigen::Vector3d> &positions,
    std::vector<Eigen::Vector3d> &velocity,
    std::vector<Eigen::Vector3d> &forces) {
  for (int i = 0; i < positions.size(); i++) {

    velocity[i] = velocity[i] + dt * (gravity + forces[i]);
    positions[i] = positions[i] + dt * velocity[i];

    if (velocity[i].norm() > v_max)
      v_max = velocity[i].norm();
  }
}