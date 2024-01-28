#include "time_integration.h"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

#include "types.h"

learnSPH::timeIntegration::semiImplicitEuler::semiImplicitEuler(
    double radius, std::vector<types::object> objects, std::vector<Eigen::Vector3d> simBoundary)
{
    this->radius = radius;
    this->v_max  = 0;
    if (simBoundary.size() > 0) {
        this->simBoundary_min = simBoundary[0];
        std::cout << "t_integ sim min: " << this->simBoundary_min << "\n";
        this->simBoundary_max = simBoundary[1];
        std::cout << "t_integ sim max: " << this->simBoundary_max << "\n";
        this->boundary_checking = true;
    } else {
        this->boundary_checking = false;
    }
    if (objects.size() > 0) {
        this->boundaries        = objects;
        this->boundary_checking = true;
    } else {
        this->boundary_checking = false;
    }
}

void learnSPH::timeIntegration::semiImplicitEuler::integrationStep(
    std::vector<Eigen::Vector3d> &positions, std::vector<Eigen::Vector3d> &velocity,
    std::vector<Eigen::Vector3d> &accelerations, std::vector<bool> &deleteFlag, double dt,
    int &count_del, Eigen::Vector3d &min_fluid_reco, Eigen::Vector3d &max_fluid_reco)
{
    v_max         = 0;
    bool copyFlag = false;
    count_del     = 0;
    if (positions.size() > 0) {
        max_fluid_reco = positions[0];
        min_fluid_reco = positions[0];
    }

    for (int i = 0; i < positions.size(); i++) {
        velocity[i]  = velocity[i] + dt * accelerations[i];
        positions[i] = positions[i] + dt * velocity[i];

        if (velocity[i].squaredNorm() > v_max)
            v_max = velocity[i].squaredNorm();

        if (positions[i].x() > max_fluid_reco.x())
            max_fluid_reco.x() = positions[i].x();
        if (positions[i].y() > max_fluid_reco.y())
            max_fluid_reco.y() = positions[i].y();
        if (positions[i].z() > max_fluid_reco.z())
            max_fluid_reco.z() = positions[i].z();
        if (positions[i].x() < min_fluid_reco.x())
            min_fluid_reco.x() = positions[i].x();
        if (positions[i].y() < min_fluid_reco.y())
            min_fluid_reco.y() = positions[i].y();
        if (positions[i].z() < min_fluid_reco.z())
            min_fluid_reco.z() = positions[i].z();

        // mark out of bound elements for deletion if activated

        if (this->boundary_checking) {
            if (positions[i].x() > simBoundary_max.x() || positions[i].x() < simBoundary_min.x() ||
                positions[i].y() > simBoundary_max.y() || positions[i].y() < simBoundary_min.y() ||
                positions[i].z() > simBoundary_max.z() || positions[i].z() < simBoundary_min.z()) {
                deleteFlag[i] = true;
                copyFlag      = true;
                count_del++;
            }

            for (int k = 0; k < boundaries.size(); k++) {
                if (positions[i].x() < boundaries[k].max.x() &&
                    positions[i].x() > boundaries[k].min.x() &&
                    positions[i].y() < boundaries[k].max.y() &&
                    positions[i].y() > boundaries[k].min.y() &&
                    positions[i].z() < boundaries[k].max.z() &&
                    positions[i].z() > boundaries[k].min.z() && !(boundaries[k].noCheck)) {
                    std::cout << "What happened here?\n";
                    deleteFlag[i] = true;
                    copyFlag      = true;
                    count_del++;
                }
            }
        }
    }
}
