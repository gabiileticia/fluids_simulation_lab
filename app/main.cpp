#include <stdlib.h>     // rand
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>    // std::max

#include <Eigen/Dense>

#include "../learnSPH/io.h"
#include "../learnSPH/kernel.h"
#include "../learnSPH/sampling.h"
#include "../learnSPH/geometry.h"

int main()
{
	std::cout << "Welcome to the learnSPH framework!!" << std::endl;
	std::cout << "Generating a sample scene...";

	// Initializing 
	double particle_radius = 0.03;
	const double particle_diameter = 2.0 * particle_radius;
	const double fluid_sampling_distance = particle_diameter;
	const double boundary_sampling_distance = 0.8 * particle_diameter;
	// smoothing lenght
	const double h = 1.2 * particle_diameter;
	// compact support
	const double beta = 2.0 * h;

	learnSPH::kernel::CubicSplineKernel cubicSpline(h);

	const double fluid_density = 1000.0;
	const double fluid_volume = 1.0;
	const double fluid_mass = fluid_volume/fluid_density;


	// Load simulation geometry

	std::vector<Eigen::Vector3d> boundary_particles;
	learnSPH::geometry::load_n_sample_boundary(boundary_particles
											   , "./res/box_density_test.obj"
											   , boundary_sampling_distance);


	// // Load a obj surface mesh
	// const std::vector<learnSPH::TriMesh> meshes = learnSPH::read_tri_meshes_from_obj("../res/box.obj");
	// const learnSPH::TriMesh& box = meshes[0];

	// // Sample the mesh with particles
	// const double sampling_distance = 0.05;
	// std::vector<Eigen::Vector3d> particles;
	// learnSPH::sampling::triangle_mesh(particles, box.vertices, box.triangles, sampling_distance);

	// // Initialize data vectors
	// std::vector<double> particles_scalar_data(particles.size());
	// std::vector<Eigen::Vector3d> particles_vector_data(particles.size());

	// // Scalar data will be the particle_id
	// for (int i = 0; i < (int)particles.size(); i++) {
	// 	particles_scalar_data[i] = i;
	// }
	
	// // Simulation loop
	// for (int time_step = 0; time_step < 100; time_step++) {

	// 	for (int particle_i = 0; particle_i < (int)particles.size(); particle_i++) {
	// 		// Move particles a bit down in the Z direction
	// 		particles[particle_i][2] -= 0.025;

	// 		// Clamp the Z coord to the floor
	// 		particles[particle_i][2] = std::max(particles[particle_i][2], -1.0);
			
	// 		// Vector data is going to be the position
	// 		particles_vector_data[particle_i] = particles[particle_i] - Eigen::Vector3d(-1, -1, -1);
	// 	}

	// 	// Save output
	// 	const std::string filename = "../res/example_" + std::to_string(time_step) + ".vtk";
	// 	learnSPH::write_particles_to_vtk(filename, particles, particles_scalar_data, particles_vector_data);
	// }

	// std::cout << "completed!" << std::endl;
	// std::cout << "The scene files have been saved in the folder `<build_folder>/res`. \nYou can visualize them with Blender by using the Blender Sequence Loaded addon." << std::endl;

	return 0;
}