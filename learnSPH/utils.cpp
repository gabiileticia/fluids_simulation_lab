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

void learnSPH::utils::deleteOutOfBounds(std::vector<Eigen::Vector3d> &positions,
                                        std::vector<Eigen::Vector3d> &velocity,
                                        std::vector<Eigen::Vector3d> &accelerations,
                                        std::vector<double> &densities,
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

Eigen::Vector3d learnSPH::utils::finiteDifference(learnSPH::types::ImplicitSurface foo,
                                                  Eigen::Vector3d x1, Eigen::Vector3d x2,
                                                  double tolerance)
{
    // Unit vectors
    Eigen::Vector3d ex = Eigen::Vector3d(1, 0, 0);
    Eigen::Vector3d ey = Eigen::Vector3d(0, 1, 0);
    Eigen::Vector3d ez = Eigen::Vector3d(0, 0, 1);

    Eigen::Vector3d fin_diff;

    fin_diff[0] = foo(x1 - x2 + tolerance * ex) - foo(x1 - x2 - tolerance * ex);
    fin_diff[1] = foo(x1 - x2 + tolerance * ey) - foo(x1 - x2 - tolerance * ey);
    fin_diff[2] = foo(x1 - x2 + tolerance * ez) - foo(x1 - x2 - tolerance * ez);

    fin_diff = fin_diff / (2.0 * tolerance);

    return fin_diff;
}