#pragma once
#include <array>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "../extern/Eigen/Eigen/Dense"
#include "types.h"

namespace learnSPH
{
namespace sampling
{
void fluid_box(std::vector<Eigen::Vector3d> &particles, const Eigen::Vector3d &bottom,
               const Eigen::Vector3d &top, const double sampling_distance);
void triangle_mesh(std::vector<Eigen::Vector3d> &particles,
                   const std::vector<Eigen::Vector3d> &vertices,
                   const std::vector<std::array<int, 3>> &triangles,
                   const double sampling_distance);

std::vector<std::array<int, 2>> _find_edges(const std::vector<std::array<int, 3>> &triangles);
void _sample_triangle(std::vector<Eigen::Vector3d> &particles, const Eigen::Vector3d &a,
                      const Eigen::Vector3d &b, const Eigen::Vector3d &c,
                      const double sampling_distance);
bool _point_on_triangle(const Eigen::Vector3d &a, const Eigen::Vector3d &b,
                        const Eigen::Vector3d &c, const Eigen::Vector3d &n,
                        const Eigen::Vector3d &p);
void fluid_spheres(std::vector<Eigen::Vector3d> &particles,
                   std::vector<Eigen::Vector3d> &velocities,
                   const std::vector<learnSPH::types::fluid_sphere> spheres,
                   const double particle_radius,
                   const std::vector<Eigen::Vector3d> fluid_velocities);
} // namespace sampling
} // namespace learnSPH