#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <random>


#include "catch.hpp"

#include "../learnSPH/kernel.h"
#include "../learnSPH/neighborhood_search.h"
#include "../learnSPH/densities.h"
#include "../learnSPH/io.h"
#include "../learnSPH/sampling.h"
#include "../extern/CompactNSearch/include/CompactNSearch/CompactNSearch"


#include <chrono>


struct DensitiesCompute
{
    bool compute_density()
    {
        return true;
    }
};


TEST_CASE( "Tests for density computation", "[density]" )
{
    DensitiesCompute densities_struct;

    double particle_radius = 0.025;
	const double particle_diameter = 2.0 * particle_radius;
	const double fluid_sampling_distance = particle_diameter;
	const double boundary_sampling_distance = 0.8 * particle_diameter;
	// smoothing lenght
	const double h = 1.2 * particle_diameter;
	// compact support
	const double beta = 2.0 * h;

	const double fluid_density = 1000.0;
	const double fluid_volume = 800000.0;
	const double fluid_mass = fluid_volume/fluid_density;

	// Load boundary surface mesh
	const std::vector<learnSPH::TriMesh> boundary_meshes = learnSPH::read_tri_meshes_from_obj("./res/box_density_test.obj");
	const learnSPH::TriMesh& boundary = boundary_meshes[0];

	// Sample the mesh with particles
	std::vector<Eigen::Vector3d> boundary_particles;
	learnSPH::sampling::triangle_mesh(boundary_particles
									, boundary.vertices
									, boundary.triangles
									, boundary_sampling_distance);

	std::cout << "Number of boundary particles" << std::endl;
	std::cout << boundary_particles.size() << std::endl;

	// Sample the fluid with particles
	std::vector<Eigen::Vector3d> particles;
	learnSPH::sampling::fluid_box(particles
								, Eigen::Vector3d(-0.8, -0.8, -1.0)
								, Eigen::Vector3d(0.8, 0.8, -0.5)
								, fluid_sampling_distance);
	
	std::cout << "Number of fluid particles" << std::endl;
	std::cout << particles.size() << std::endl;

	double fluid_particle_mass = fluid_mass/particles.size();

	// Neighborhood search
	CompactNSearch::NeighborhoodSearch nsearch(beta);
	unsigned int point_set_id_boundary = nsearch.add_point_set(boundary_particles.front().data()
															 , boundary_particles.size());

	unsigned int point_set_id_fluid = nsearch.add_point_set(particles.front().data()
														  , particles.size());
														
	nsearch.set_active(point_set_id_boundary, point_set_id_fluid, false);
	nsearch.set_active(point_set_id_fluid, point_set_id_fluid, true);
	nsearch.set_active(point_set_id_fluid, point_set_id_boundary, true);
	nsearch.find_neighbors();

	CompactNSearch::PointSet const& ps_boundary = nsearch.point_set(point_set_id_boundary);
	CompactNSearch::PointSet const& ps_fluid = nsearch.point_set(point_set_id_fluid);

	// Compute boundary masses
	std::vector<double> boundary_particles_masses(boundary_particles.size());
	learnSPH::densities::compute_boundary_masses(boundary_particles_masses, boundary_particles, point_set_id_boundary, ps_boundary, fluid_density, h);

    // Compute fluid particles densities
	std::vector<double> particles_densities(particles.size());
	learnSPH::densities::compute_fluid_density(particles_densities
											, particles
											, boundary_particles
											, boundary_particles_masses
											, point_set_id_fluid
											, ps_fluid
											, point_set_id_boundary
											, ps_boundary
											, fluid_particle_mass
											, h);

	std::vector<double> boundary_particles_densities(boundary_particles.size(), fluid_density);
	boundary_particles.insert(boundary_particles.end(), particles.begin(), particles.end());
	boundary_particles_densities.insert(boundary_particles_densities.end(), particles_densities.begin(), particles_densities.end());
	// Save output
	const std::string filename = "./res/assign2.vtk";
	learnSPH::write_particles_to_vtk(filename, boundary_particles, boundary_particles_densities);

	// const std::string filename = "./res/assign2.vtk";
	// learnSPH::write_particles_to_vtk(filename, boundary_particles, boundary_particles_masses);

    SECTION("--") {
        double sum_mass = 0.0;
        double max = -10;
        double min = 10;
        for(int i=0; i < boundary_particles_masses.size(); ++i){
            sum_mass += boundary_particles_masses[i];
            if (boundary_particles_masses[i] > max){
                max = boundary_particles_masses[i];
            }
            if (boundary_particles_masses[i] < min){
                min = boundary_particles_masses[i];
            }
        }
        std::cout << sum_mass/boundary_particles_masses.size() << std::endl;
        std::cout << sum_mass << std::endl;
        std::cout << max << std::endl;
        std::cout << min << std::endl;
    }
}
