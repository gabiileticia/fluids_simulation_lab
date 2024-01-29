#include "utils.h"
#include "kernel.h"

#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdlib.h> // rand
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

int learnSPH::utils::deleteOutOfBounds(std::vector<Eigen::Vector3d> &positions,
                                       std::vector<Eigen::Vector3d> &velocity,
                                       std::vector<Eigen::Vector3d> &accelerations,
                                       std::vector<double> &densities,
                                       std::vector<double> &pressure, std::vector<bool> &deleteFlag,
                                       int &count_del)
{

    // std::cout << "Deleting " << count_del << " elements from particle vectors." << std::endl;
    int counter = 0;
    for (int i = 0; i < positions.size(); i++) {
        // skip marked for deletion element and don't increase conter
        if (deleteFlag[i]) {
            continue;
        }
        // copy only if element has to be shifted
        if (counter != i) {
            positions[counter]     = positions[i];
            velocity[counter]      = velocity[i];
            accelerations[counter] = accelerations[i];
            densities[counter]     = densities[i];
            pressure[counter]      = pressure[i];
        }
        counter++;
    }
    positions.resize(counter);
    velocity.resize(counter);
    accelerations.resize(counter);
    densities.resize(counter);
    pressure.resize(counter);
    deleteFlag.resize(counter);
    std::fill(deleteFlag.begin(), deleteFlag.end(), false);

    return count_del;
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
    oss << std::put_time(timeInfo, "%m-%d_%H-%M-%S");
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

std::ostringstream learnSPH::utils::updateProgressBar(int currentStep, int maxSteps,
                                                      const int barWidth)
{
    std::ostringstream output_msg;

    static std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();

    float progress        = static_cast<float>(currentStep) / maxSteps;
    int progressBarLength = static_cast<int>(progress * barWidth);

    output_msg << "Frame: " << currentStep << "/" << maxSteps;
    output_msg << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < progressBarLength) {
            output_msg << "#";
        } else {
            output_msg << " ";
        }
    }
    output_msg << "] " << int(progress * 100.0) << "% | ";

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

    output_msg << "Elapsed Time: " << std::setfill('0') << std::setw(2) << elapsedHours << ":"
               << std::setfill('0') << std::setw(2) << elapsedMins << ":" << std::setfill('0')
               << std::setw(2) << elapsedSecs;
    output_msg << " | Remaining Time: " << std::setfill('0') << std::setw(2) << remainingHours
               << ":" << std::setfill('0') << std::setw(2) << remainingMins << ":"
               << std::setfill('0') << std::setw(2) << remainingSecs;

    return output_msg;
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

void learnSPH::utils::logMessage(const std::string &message, const std::string &filename)
{
    std::ofstream logfile;
    logfile.open(filename, std::ios_base::app); // Open file in append mode

    if (logfile.is_open()) {
        logfile << message << std::endl; // Write message to the file
        logfile.close();                 // Close the file
    } else {
        std::cerr << "Error opening the file." << std::endl;
    }
}

double learnSPH::utils::particle_mass(const double fluid_rest_density,
                                      const double sampling_distance)
{
    double pack_density, particle_volume, fluid_mass, sample_volume;

    // sample volume is 1 cubic meter
    sample_volume   = 1.0;
    pack_density    = learnSPH::kernel::PI * std::sqrt(2) / 6;
    particle_volume = (1.0 / 6.0) * learnSPH::kernel::PI * sampling_distance * sampling_distance *
                      sampling_distance;
    fluid_mass = 1.0 * fluid_rest_density;

    return fluid_mass / (sample_volume * pack_density / particle_volume);
}

void learnSPH::utils::create_emitter_shield(const Eigen::Matrix3d &rotationMatrix,
                                            const Eigen::Vector3d emitOrigin,
                                            const double emitRadius,
                                            std::vector<Eigen::Vector3d> &boundaryParticles,
                                            const double particleRadius,
                                            unsigned int &point_set_id_boundary,
                                            CompactNSearch::NeighborhoodSearch &nsearch)
{
    double shieldRadius, squaredShieldRadius, corner, x, y, z, particleDiameter, z_origin;
    int l, new_part_counter, old_boundary_size;

    std::vector<Eigen::Vector3d> shield_particles;

    particleDiameter      = particleRadius * 2;
    shieldRadius          = emitRadius + 1.5 * particleDiameter;
    new_part_counter      = 0;
    l                     = (int)(shieldRadius / particleRadius);
    x                     = 0;
    y                     = 0;
    z                     = 0;
    z_origin              = -1.5 * particleDiameter;
    corner                = -l * particleRadius;
    squaredShieldRadius   = l * particleRadius * l * particleRadius;
    const double sqrt3    = std::sqrt(3);
    const double onethird = 1. / 3.;
    const double twothird = 2. / 3.;
    const double sqrt6    = std::sqrt(6);
    const double epsilon  = 5e-4;

    for (int k = 0; k < 4; k++) {
        z = z_origin + particleRadius * k * twothird * sqrt6;
        for (int i = 0; i < 2 * l; i++) {
            for (int j = 0; j < 2 * l; j++) {
                x = corner + particleRadius * (2 * i + ((j + k) % 2));
                y = corner + particleRadius * (sqrt3 * (j + onethird * (k % 2)));
                if ((k == 0) && ((x * x + y * y - squaredShieldRadius) < 0)) {
                    shield_particles.push_back(rotationMatrix * Eigen::Vector3d({x, y, z}) +
                                               emitOrigin);
                    new_part_counter++;
                } else if ((k != 0) &&
                           (std::abs((x * x + y * y - squaredShieldRadius)) < epsilon)) {
                    shield_particles.push_back(rotationMatrix * Eigen::Vector3d({x, y, z}) +
                                               emitOrigin);
                    new_part_counter++;
                }
            }
        }
    }
    old_boundary_size = boundaryParticles.size();

    boundaryParticles.resize(old_boundary_size + new_part_counter);

    nsearch.resize_point_set(point_set_id_boundary, boundaryParticles.front().data(),
                             boundaryParticles.size());

    for (int i = 0; i < new_part_counter; i++) {
        boundaryParticles[old_boundary_size + i] = shield_particles[i];
    }
}

void learnSPH::utils::zeroCheck(std::vector<Eigen::Vector3d> &particles, std::string msg,
                                double epsilon)
{
    for (int i = 0; i < particles.size(); i++) {
        if ((particles[i].norm() - 0.0) < epsilon) {
            std::cout << msg;
            exit(-1);
        }
    }
}

void learnSPH::utils::checkBoundaryLifetimes(
    std::vector<Eigen::Vector3d> &boundary_particles, unsigned int &ps_id_boundary,
    CompactNSearch::NeighborhoodSearch &nsearch,
    std::vector<std::array<int, 2>> &boundary_particle_info,
    std::vector<learnSPH::types::object> &object_info, double currentTime)
{
    for (int i = 0; i < object_info.size(); i++) {
        if ((object_info[i].lifetime <= currentTime) && (object_info[i].lifetime != 0)) {
            boundary_particles.erase(boundary_particles.begin() + boundary_particle_info[i][0],
                                     boundary_particles.begin() + boundary_particle_info[i][1]);
                                     
            boundary_particle_info.erase(boundary_particle_info.begin() + i);
            object_info.erase(object_info.begin() + i);
        }
    }
    nsearch.resize_point_set(ps_id_boundary, boundary_particles.front().data(),
                             boundary_particles.size());

    nsearch.find_neighbors();
}