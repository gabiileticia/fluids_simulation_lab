#include "utils.h"

#include <array>
#include <cerrno>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <fstream>
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

int learnSPH::utils::cubeVertex2VertexIndex(uint cellIdx, uint vertexIndex, uint nx, uint ny,
                                            uint nz)
{
    uint i, j, k;
    k = cellIdx % nz;
    j = (cellIdx / nz) % ny;
    i = cellIdx / (ny * nz);

    switch (vertexIndex) {
    case 0:
        break;
    case 1:
        i += 1;
        break;
    case 2:
        i += 1;
        j += 1;
        break;
    case 3:
        j += 1;
        break;
    case 4:
        k += 1;
        break;
    case 5:
        i += 1;
        k += 1;
        break;
    case 6:
        i += 1;
        j += 1;
        k += 1;
        break;
    case 7:
        j += 1;
        k += 1;
        break;
    default:
        std::cout << "vertex2VertexIdx: Mismatch of lookup index: " << vertexIndex << "\n";
        exit(-1);
    }
    if (0 <= i < nx && 0 <= j < ny && 0 <= k < nz) {
        return i * ny * nz + j * nz + k;
    } else {
        return -1;
    }
}

/**
 * Neighbor mapping from center vertex v:
 *
 *          2   1
 *          |  /
 *          | /
 *          |/
 * 3 ______ v ______ 0
 *         /|
 *        / |
 *       /  |
 *      4   5
 */

int learnSPH::utils::vertexSixNeighbors(uint vertexIndex, int neighbor, uint nx, uint ny, uint nz)
{
    uint i, j, k;
    k = vertexIndex % nz;
    j = (vertexIndex / nz) % ny;
    i = vertexIndex / (ny * nz);

    switch (neighbor) {
    case 0:
        i += 1;
        break;
    case 1:
        j += 1;
        break;
    case 2:
        k += 1;
        break;
    case 3:
        i -= 1;
        break;
    case 4:
        j -= 1;
        break;
    case 5:
        k -= 1;
        break;

    default:
        std::cout << "sixNeighbors: Mismatch of lookup index: " << vertexIndex << "\n";
        exit(-1);
    }
    if (0 <= i < nx && 0 <= j < ny && 0 <= k < nz) {
        return i * ny * nz + j * nz + k;
    } else {
        return -1;
    }
}

int64_t learnSPH::utils::vertex8NeighborCells(uint cellIndex, int neighbor, uint nx, uint ny,
                                              uint nz)
{
    uint i, j, k;
    k = cellIndex % nz;
    j = (cellIndex / nz) % ny;
    i = cellIndex / (ny * nz);

    switch (neighbor) {
    case 0:
        i -= 1;
        j -= 1;
        k -= 1;
        break;
    case 1:
        j -= 1;
        k -= 1;
        break;
    case 2:
        k -= 1;
        break;
    case 3:
        i -= 1;
        k -= 1;
        break;
    case 4:
        i -= 1;
        j -= 1;
        break;
    case 5:
        j -= 1;
        break;
    case 6:
        break;
    case 7:
        i -= 1;
        break;
    default:
        exit(-1);
    }

    if (0 <= i < nx && 0 <= j < ny && 0 <= k < nz) {
        return i * ny * nz + j * nz + k;
    } else {
        return -1;
    }
}

std::array<int, 4> learnSPH::utils::celladjByEdge(int edge)
{
    switch (edge) {
    case 0:
        return {1, 2, 5, 6};
    case 1:
        return {2, 3, 6, 7};
    case 2:
        return {4, 5, 6, 7};
    case 3:
        return {0, 3, 4, 7};
    case 4:
        return {0, 1, 4, 5};
    case 5:
        return {0, 1, 2, 3};
    default:
        exit(-1);
    }
}

Eigen::Vector3d learnSPH::utils::index2coord(uint vertexIndex, double cellwidth, uint nx, uint ny,
                                             uint nz, Eigen::Vector3d origin)
{
    int x                  = vertexIndex / (ny * nz);
    int y                  = (vertexIndex / nz) % ny;
    int z                  = vertexIndex % nz;
    Eigen::Vector3d coords = Eigen::Vector3d(x * cellwidth, y * cellwidth, z * cellwidth) + origin;
    return coords;
}

void learnSPH::utils::logMessage(const std::string& message, const std::string& filename) {
    std::ofstream logfile(filename);
    logfile.open(filename, std::ios_base::app); // Open file in append mode

    if (logfile.is_open()) {
        logfile << message << std::endl; // Write message to the file
        logfile.close(); // Close the file
    } else {
        std::cerr << "Error opening the file." << std::endl;
    }
}