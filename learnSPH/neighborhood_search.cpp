#include "neighborhood_search.h"
#include <iostream>

std::vector<std::string> learnSPH::neighborhood_search::brute_force_neighbor(const double h, std::vector<Eigen::Vector3d> particles, int particle_index, double beta)
{
	// std::cout << "current particle" << std::endl;
	// std::cout << particles[particle_index] << std::endl;

	std::vector<std::string> particle_neighbors(particles.size());

	// std::cout << "neighbors" << std::endl;
	for (int i = 0; i < (int)particles.size(); i++) {
		if ((particles[i] - particles[particle_index]).norm() <= beta * h){
			particle_neighbors[particle_index] = particle_neighbors[particle_index] + ';' + std::to_string(i) ;
		}
	}
	// std::cout << particle_neighbors[particle_index] << std::endl;
	return particle_neighbors;
}

void learnSPH::neighborhood_search::brute_force_neighbors(const double h, std::vector<Eigen::Vector3d> particles, double beta){
	for (int i = 0; i < (int)particles.size(); i++) {
		brute_force_neighbor(h, particles, i, beta);
	}
}
