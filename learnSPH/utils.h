#ifndef UTILS
#define UTILS

#include <Eigen/Dense>
#include <array>
#include <chrono>
#include <cstdint>
#include <vector>

#include "types.h"

namespace learnSPH
{
    namespace utils
    {
        int deleteOutOfBounds(std::vector<Eigen::Vector3d> &positions,
                            std::vector<Eigen::Vector3d> &velocity,
                            std::vector<Eigen::Vector3d> &accelerations,
                            std::vector<double> &densities,
                            std::vector<double> &pressure, std::vector<bool> &deleteFlat,
                            int &count_del);

        void create_simulation_folder(const std::string assign_number, std::string &timestamp);

        void updateProgressBar(int currentStep, int maxSteps, const int barWidth);

        Eigen::Vector3d implicitVertexNormal(learnSPH::types::ImplicitSurface foo, Eigen::Vector3d vertex,
                                            double epsilon, void *fooArgs);

        int cubeVertex2VertexIndex(uint cellIdx, uint vertexIndex, uint nx, uint ny, uint nz);

        int vertexSixNeighbors(uint vertexIndex, int neighbor, uint nx, uint ny, uint nz);

        int64_t vertex8NeighborCells(uint cellIndex, int neighbor, uint nx, uint ny, uint nz);

        Eigen::Vector3d index2coord(uint vertexIndex, double cellwidth, uint nx , uint ny, uint nz, Eigen::Vector3d origin);

        std::array<int, 4> celladjByEdge(int edge);

        void logMessage(const std::string& message, const std::string& filename);
    } // namespace utils
} // namespace learnSPH

#endif