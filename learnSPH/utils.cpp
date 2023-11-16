#include "utils.h"
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <stdlib.h> // rand
#include <iostream>
#include <chrono>
#include <sys/stat.h>


void learnSPH::utils::deleteOutOfBounds(std::vector<Eigen::Vector3d> &positions,
                                        std::vector<Eigen::Vector3d> &velocity,
                                        std::vector<Eigen::Vector3d> &accelerations,
                                        std::vector<double> &densities,
                                        std::vector<double> &pressure,
                                        std::vector<bool> &deleteFlag,
                                        int &count_del)
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


void learnSPH::utils::create_simulation_folder(const std::string assign_number, std::string &timestamp)
{
    std::time_t result = std::time(nullptr);
    std::stringstream strm;
    strm << result;

    timestamp = strm.str();

    // Create folder
    // doesn't work
    std::string stringpath = "./res/" + assign_number + "/" + timestamp + "/";
    int status = mkdir(stringpath.c_str(),0777);
    if (status == -1){
        std::cerr << "Error: " << strerror(errno) << "\n";
        exit(errno);
    }
    else 
        return;
    
    // Create file with setup


}

