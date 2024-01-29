#ifndef UTILS
#define UTILS

#include <array>
#include <chrono>
#include <cstdint>
#include <sstream>
#include <vector>

#include "../CompactNSearch/include/CompactNSearch/CompactNSearch.h"
#include "../extern/Eigen/Eigen/Dense"
#include "types.h"

namespace learnSPH
{
namespace utils
{
int deleteOutOfBounds(std::vector<Eigen::Vector3d> &positions,
                      std::vector<Eigen::Vector3d> &velocity,
                      std::vector<Eigen::Vector3d> &accelerations, std::vector<double> &densities,
                      std::vector<double> &pressure, std::vector<bool> &deleteFlat, int &count_del);

void create_simulation_folder(const std::string assign_number, std::string &timestamp);

std::ostringstream updateProgressBar(int currentStep, int maxSteps, const int barWidth);

Eigen::Vector3d implicitVertexNormal(learnSPH::types::ImplicitSurface foo, Eigen::Vector3d vertex,
                                     double epsilon, void *fooArgs);

int cubeVertex2VertexIndex(uint cellIdx, uint vertexIndex, uint nx, uint ny, uint nz);

int vertexSixNeighbors(uint vertexIndex, int neighbor, uint nx, uint ny, uint nz);

int64_t vertex8NeighborCells(uint cellIndex, int neighbor, uint nx, uint ny, uint nz);

Eigen::Vector3d index2coord(uint vertexIndex, double cellwidth, uint nx, uint ny, uint nz,
                            Eigen::Vector3d origin);

std::array<int, 4> celladjByEdge(int edge);

void logMessage(const std::string &message, const std::string &filename);
double particle_mass(const double fluid_rest_density, const double sampling_distance);

void create_emitter_shield(const Eigen::Matrix3d &rotationMatrix, const Eigen::Vector3d emitOrigin,
                           const double emitRadius, std::vector<Eigen::Vector3d> &boundaryParticles,
                           const double particlesRadius, unsigned int &point_set_id_boundary,
                           CompactNSearch::NeighborhoodSearch &nsearch);
void zeroCheck(std::vector<Eigen::Vector3d> &particles, std::string msg, double epsilon);

void checkBoundaryLifetimes(std::vector<Eigen::Vector3d> &boundary_particles,
                            unsigned int &ps_id_boundary,
                            CompactNSearch::NeighborhoodSearch &nsearch,
                            std::vector<std::array<int, 2>> &boundary_info,
                            std::vector<learnSPH::types::object> &object_info,
                            double currentTime);
} // namespace utils
} // namespace learnSPH

#endif