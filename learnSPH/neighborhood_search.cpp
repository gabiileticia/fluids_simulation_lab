#include "neighborhood_search.h"
#include <iostream>

std::vector<std::string> learnSPH::neighborhood_search::brute_force_neighbor(const double h, std::vector<Eigen::Vector3d> particles, int particle_index,double beta) {

  std::vector<std::string> particle_neighbors(particles.size());

  for (int i = 0; i < (int)particles.size(); i++) {
    if ((particles[i] - particles[particle_index]).norm() <= beta * h) {
      particle_neighbors[particle_index] =
          particle_neighbors[particle_index] + ';' + std::to_string(i);
    }
  }

  return particle_neighbors;
}

void learnSPH::neighborhood_search::brute_force_neighbors(const double h, std::vector<Eigen::Vector3d> particles, double beta) {
  for (int i = 0; i < (int)particles.size(); i++) {
    brute_force_neighbor(h, particles, i, beta);
  }
}

std::vector<std::vector<int>> learnSPH::neighborhood_search::brute_force_search(const double h, std::vector<Eigen::Vector3d> particles, double beta) {

		std::vector<std::vector<int>> neighbors(particles.size());
		double distance;
		
		for (int i = 0; i < particles.size(); ++i){
			for (int j = 0; j < particles.size(); j++){
				if (i != j){
					distance = (particles[i] - particles[j]).norm();
					if(distance <= beta * h){
						neighbors[i].push_back(j);
					}
				}
			}
		}
		return neighbors;
}