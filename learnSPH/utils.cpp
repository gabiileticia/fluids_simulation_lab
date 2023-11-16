#include "utils.h"

#include <cerrno>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iostream>
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

void learnSPH::utils::updateProgressBar(int &currentStep, int &maxSteps) {
    static std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();

    const int barWidth = 50;
    float progress = static_cast<float>(currentStep) / maxSteps;
    int progressBarLength = static_cast<int>(progress * barWidth);

    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < progressBarLength) {
            std::cout << "=";
        } else {
            std::cout << " ";
        }
    }
    std::cout << "] " << int(progress * 100.0) << "%\r";
    std::cout.flush();

    // Calculate elapsed time
    std::chrono::steady_clock::time_point currentTime = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsedSeconds = currentTime - startTime;
    
    // Calculate remaining time
    double remainingSeconds = elapsedSeconds.count() / progress - elapsedSeconds.count();

    std::cout << "Elapsed Time: " << elapsedSeconds.count() << " seconds";
    std::cout << " | Remaining Time: " << remainingSeconds << " seconds    ";
}