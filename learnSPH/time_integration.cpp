#include "time_integration.h"

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <vector>

#include "types.h"

learnSPH::timeIntegration::semiImplicitEuler::semiImplicitEuler(double dt,
                                                                double radius,
                                                                double v_max) {
  assert(dt <=
         0.5 * radius / std::pow(std::pow(v_max, 2) + std::pow(v_max, 2), 0.5));
  dt = dt;
  vec_dt = Eigen::Vector3d(dt, dt, dt);
  radius = radius;
  v_max = v_max;
}

void learnSPH::timeIntegration::semiImplicitEuler::integrationStep(
    std::vector<learnSPH::types::Particle> &particles) {
  for (int i = 0; i < particles.size(); i++) {
    particles[i].velocity = particles[i].velocity + vec_dt * gravity;
    particles[i].pos = particles[i].pos + dt * particles[i].velocity;
  }
}