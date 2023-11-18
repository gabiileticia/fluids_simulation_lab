#include "time_integration.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

#include "types.h"

learnSPH::timeIntegration::semiImplicitEuler::semiImplicitEuler(double radius)
{
    this->radius            = radius;
    this->v_max             = 0;
    this->min_boundary      = {0, 0, 0};
    this->max_boundary      = {0, 0, 0};
    this->boundary_checking = false;
}
learnSPH::timeIntegration::semiImplicitEuler::semiImplicitEuler(double radius,
                                                                Eigen::Vector3d max_boundary,
                                                                Eigen::Vector3d min_boundary)
{
    this->radius            = radius;
    this->v_max             = 0;
    this->min_boundary      = min_boundary;
    this->max_boundary      = max_boundary;
    this->boundary_checking = true;
}

void learnSPH::timeIntegration::semiImplicitEuler::integrationStep(
    std::vector<Eigen::Vector3d> &positions, std::vector<Eigen::Vector3d> &velocity,
    std::vector<Eigen::Vector3d> &accelerations, std::vector<bool> &deleteFlag, double dt, int &count_del)
{
    v_max         = 0;
    bool copyFlag = false;
    count_del = 0;

    for (int i = 0; i < positions.size(); i++) {
        velocity[i]  = velocity[i] + dt * accelerations[i];
        positions[i] = positions[i] + dt * velocity[i];

        if (velocity[i].squaredNorm() > v_max)
            v_max = velocity[i].squaredNorm();

        // mark out of bound elements for deletion if activated
        if (this->boundary_checking &&
            (positions[i].x() > max_boundary.x() || positions[i].x() < min_boundary.x() ||
             positions[i].y() > max_boundary.y() || positions[i].y() < min_boundary.y() ||
             positions[i].z() > max_boundary.z() || positions[i].z() < min_boundary.z())) {
            deleteFlag[i] = true;
            copyFlag      = true;
            count_del++;
        }
    }
}