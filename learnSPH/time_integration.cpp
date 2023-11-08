#include "time_integration.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

#include "types.h"

learnSPH::timeIntegration::semiImplicitEuler::semiImplicitEuler(double radius) {
  this->radius = radius;
  this->v_max = 0;
  this->min_boundary = {0,0,0};
  this->max_boundary = {0,0,0};
  this->boundary_checking = false;
}
learnSPH::timeIntegration::semiImplicitEuler::semiImplicitEuler(
    double radius, Eigen::Vector3d max_boundary, Eigen::Vector3d min_boundary) {
  this->radius = radius;
  this->v_max = 0;
  this->min_boundary = min_boundary;
  this->max_boundary = max_boundary;
  this->boundary_checking = true;
}

void learnSPH::timeIntegration::semiImplicitEuler::integrationStep(
    std::vector<Eigen::Vector3d> &positions,
    std::vector<Eigen::Vector3d> &velocity,
    std::vector<Eigen::Vector3d> &accelerations, double dt) {
  v_max = 0;
  bool copyFlag = false;
  int count_del = 0;
  std::vector<bool> deleteFlat(positions.size(), false);

  for (int i = 0; i < positions.size(); i++) {
    velocity[i] = velocity[i] + dt * (accelerations[i]);
    positions[i] = positions[i] + dt * velocity[i];

    if (velocity[i].norm() > v_max)
      v_max = velocity[i].norm();

    if (this->boundary_checking && (positions[i].x() > max_boundary.x() ||
                                    positions[i].x() < min_boundary.x() ||
                                    positions[i].y() > max_boundary.y() ||
                                    positions[i].y() < min_boundary.y() ||
                                    positions[i].z() > max_boundary.z() ||
                                    positions[i].z() < min_boundary.z())) {
      deleteFlat[i] = true;
      copyFlag = true;
      count_del++;
    }
  }

  if (copyFlag && boundary_checking) {
    std::cout << "Deleting " << count_del << " elements from particle vectors."
              << std::endl;
    int counter = 0;
    for (int i = 0; i < positions.size(); i++) {
      // skip marked for deletion element and don't increase conter
      if (deleteFlat[i]) {
        continue;
      }
      // copy only if element has to be shifted
      if(counter != i){
        positions[counter] = positions[i];
        velocity[counter] = velocity[i];
        accelerations[counter] = accelerations[i];
      }
      counter++;
    }
    positions.resize(counter);
    velocity.resize(counter);
    accelerations.resize(counter);

    std::cout << "Done. New number of particles is: " << positions.size()
            << std::endl;
    fflush(stdout);
  }
  
}