#include "utils.h"
#include <iostream>

void learnSPH::utils::deleteOutOfBounds(std::vector<Eigen::Vector3d> &positions,
                                        std::vector<Eigen::Vector3d> &velocity,
                                        std::vector<Eigen::Vector3d> &accelerations,
                                        std::vector<double> &densities,
                                        std::vector<double> &pressure,
                                        std::vector<bool> &deleteFlag)
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

    std::cout << "Done. New number of particles is: " << positions.size() << std::endl;
    fflush(stdout);
}
}