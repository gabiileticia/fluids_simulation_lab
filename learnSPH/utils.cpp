#include "utils.h"

#include <array>
#include <cerrno>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdlib.h> // rand
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

void learnSPH::utils::deleteOutOfBounds(std::vector<Eigen::Vector3d> &positions,
                                        std::vector<Eigen::Vector3d> &velocity,
                                        std::vector<Eigen::Vector3d> &accelerations,
                                        std::vector<double> &densities,
                                        std::vector<double> &fluid_reco_densities,
                                        std::vector<double> &pressure,
                                        std::vector<bool> &deleteFlag, int &count_del)
{

    std::cout << "Deleting " << count_del << " elements from particle vectors." << std::endl;
    int counter = 0;
    for (int i = 0; i < positions.size(); i++) {
        // skip marked for deletion element and don't increase conter
        if (deleteFlag[i]) {
            continue;
        }
        // copy only if element has to be shifted
        if (counter != i) {
            positions[counter]            = positions[i];
            velocity[counter]             = velocity[i];
            accelerations[counter]        = accelerations[i];
            densities[counter]            = densities[i];
            fluid_reco_densities[counter] = fluid_reco_densities[i];
            pressure[counter]             = pressure[i];
        }
        counter++;
    }
    positions.resize(counter);
    velocity.resize(counter);
    accelerations.resize(counter);
    densities.resize(counter);
    fluid_reco_densities.resize(counter);
    pressure.resize(counter);
    deleteFlag.resize(counter);
    std::fill(deleteFlag.begin(), deleteFlag.end(), false);

    std::cout << "Done. New number of particles is: " << positions.size() << std::endl;
    fflush(stdout);
}

void learnSPH::utils::create_simulation_folder(const std::string assign_number,
                                               std::string &timestamp)
{
    // Get current time
    auto currentTime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    // Convert current time to a struct tm for formatting
    std::tm *timeInfo = std::localtime(&currentTime);

    // Check if conversion failed
    if (timeInfo == nullptr) {
        std::cerr << "Failed to get local time.\n";
        exit(errno);
    }

    // Format the time into HH:MM:SS
    std::ostringstream oss;
    oss << std::put_time(timeInfo, "%H_%M_%S");
    timestamp = oss.str();

    // Create folder
    // doesn't work
    std::string stringpath = "./res/" + assign_number + "/" + timestamp + "/";
    int status             = mkdir(stringpath.c_str(), 0777);

    try {
        std::filesystem::create_directories(stringpath);
    } catch (const std::exception e) {
        std::cerr << "Error: " << e.what() << "\n";
        exit(errno);
    }
}

void learnSPH::utils::updateProgressBar(int currentStep, int maxSteps, const int barWidth)
{
    static std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();

    float progress        = static_cast<float>(currentStep) / maxSteps;
    int progressBarLength = static_cast<int>(progress * barWidth);

    std::cout << "Frame: " << currentStep << "/" << maxSteps;
    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < progressBarLength) {
            std::cout << "#";
        } else {
            std::cout << " ";
        }
    }
    std::cout << "] " << int(progress * 100.0) << "% | ";

    // Calculate elapsed time
    std::chrono::steady_clock::time_point currentTime = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsedSeconds      = currentTime - startTime;

    // Calculate remaining time
    double remainingSeconds = (elapsedSeconds.count() / currentStep) * (maxSteps - currentStep);

    int elapsedHours = static_cast<int>(elapsedSeconds.count()) / 3600;
    int elapsedMins  = static_cast<int>(elapsedSeconds.count()) / 60 % 60;
    int elapsedSecs  = static_cast<int>(elapsedSeconds.count()) % 60;

    int remainingHours = static_cast<int>(remainingSeconds) / 3600;
    int remainingMins  = static_cast<int>(remainingSeconds) / 60 % 60;
    int remainingSecs  = static_cast<int>(remainingSeconds) % 60;

    std::cout << "Elapsed Time: " << std::setfill('0') << std::setw(2) << elapsedHours << ":"
              << std::setfill('0') << std::setw(2) << elapsedMins << ":" << std::setfill('0')
              << std::setw(2) << elapsedSecs;
    std::cout << " | Remaining Time: " << std::setfill('0') << std::setw(2) << remainingHours << ":"
              << std::setfill('0') << std::setw(2) << remainingMins << ":" << std::setfill('0')
              << std::setw(2) << remainingSecs << "\n";

    std::cout.flush();
}

/* Cube description:
 *         7 ________ 6
 *         /|       /|
 *       /  |     /  |
 *   4 /_______ /    |
 *    |     |  |5    |
 *    |    3|__|_____|2
 *    |    /   |    /
 *    |  /     |  /
 *    |/_______|/
 *   0          1
 */

int learnSPH::utils::cubeVertex2VertexIndex(uint cellIdx, uint vertexIndex,
                                            std::vector<Eigen::Vector3d> &gridVertices, uint nx,
                                            uint ny, uint nz)
{
    int retIdx;
    switch (vertexIndex) {
    case 0:
        return cellIdx;
    case 1:
        retIdx = cellIdx + ny * nz;
        break;
    case 2:
        retIdx = cellIdx + ny * nz + nz;
        break;
    case 3:
        retIdx = cellIdx + nz;
        break;
    case 4:
        retIdx = cellIdx + 1;
        break;
    case 5:
        retIdx = cellIdx + 1 + ny * nz;
        break;
    case 6:
        retIdx = cellIdx + 1 + ny * nz + nz;
        break;
    case 7:
        retIdx = cellIdx + 1 + nz;
        break;
    default:
        std::cout << "vertex2VertexIdx: Mismatch of lookup index: " << vertexIndex << "\n";
        exit(-1);
    }
    if (cellIdx > retIdx || retIdx > nx * ny * nz ||
        gridVertices[cellIdx].y() > gridVertices[retIdx].y() ||
        gridVertices[cellIdx].z() > gridVertices[retIdx].z())
        return -1;
    else
        return retIdx;
}

/**
 * Neighbor mapping from center vertex v:
 *
 *          5   1
 *          |  /
 *          | /
 *          |/
 * 2 ______ v ______ 0
 *         /|
 *        / |
 *       /  |
 *      3   4
 */

int learnSPH::utils::vertexSixNeighbors(uint vertexIndex, int neighbor,
                                        std::vector<Eigen::Vector3d> &gridVertices, uint nx,
                                        uint ny, uint nz)
{
    int retidx;
    int max = nx * ny * nz;
    switch (neighbor) {
    case 0:
        retidx = vertexIndex + ny * nz;
        if (vertexIndex > retidx || retidx > max)
            return -1;
        else
            return vertexIndex + ny * nz;
    case 1:
        retidx = vertexIndex + nz;
        if (vertexIndex > retidx || retidx > max ||
            gridVertices[vertexIndex].y() > gridVertices[retidx].y())
            return -1;
        else
            return vertexIndex + nz;
    case 2:
        retidx = vertexIndex - ny * nz;
        if (vertexIndex < retidx || retidx > max)
            return -1;
        else
            return vertexIndex - ny * nz;
    case 3:
        retidx = vertexIndex - nz;
        if (vertexIndex < retidx || retidx > max ||
            gridVertices[vertexIndex].y() < gridVertices[retidx].y())
            return -1;
        else
            return vertexIndex - nz;
    case 4:
        retidx = vertexIndex - 1;
        if (vertexIndex < retidx || retidx > max ||
            gridVertices[vertexIndex].z() < gridVertices[retidx].z())
            return -1;
        else
            return vertexIndex - 1;
    case 5:
        retidx = vertexIndex + 1;
        if (vertexIndex > retidx || retidx > max ||
            gridVertices[vertexIndex].z() > gridVertices[retidx].z())
            return -1;
        else
            return vertexIndex + 1;
    default:
        std::cout << "sixNeighbors: Mismatch of lookup index: " << vertexIndex << "\n";
        exit(-1);
    }
}

int64_t learnSPH::utils::vertex8NeighborCells(uint cellIndex, int neighbor, uint nx, uint ny,
                                              uint nz)
{
    int64_t retidx;
    switch (neighbor) {
    case 0:
        retidx = cellIndex - 1 - ny * nz - nz;
        break;
    case 1:
        retidx = cellIndex - 1 - nz;
        break;
    case 2:
        retidx = cellIndex - 1;
        break;
    case 3:
        retidx = cellIndex - 1 - ny * nz;
        break;
    case 4:
        retidx = cellIndex - ny * nz - nz;
        break;
    case 5:
        retidx = cellIndex - nz;
        break;
    case 6:
        retidx = cellIndex;
        break;
    case 7:
        retidx = cellIndex - ny * nz;
        break;
    default:
        exit(-1);
    }

    // if (cellIndex < retidx || retidx > nx * ny * nz)
    //     return -1;
    // else
        return retidx;
}

std::array<int, 4> learnSPH::utils::celladjByEdge(int edge)
{
    switch (edge) {
    case 0:
        return {1, 2, 5, 6};
    case 1:
        return {2, 3, 6, 7};
    case 2:
        return {0, 3, 4, 7};
    case 3:
        return {0, 1, 4, 5};
    case 4:
        return {0, 1, 2, 3};
    case 5:
        return {4, 5, 6, 7};
    default:
        exit(-1);
    }
}