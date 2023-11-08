#include <algorithm> // std::max
#include <iostream>
#include <ostream>
#include <stdlib.h> // rand
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "../extern/CompactNSearch/include/CompactNSearch/CompactNSearch"
#include "../learnSPH/acceleration.h"
#include "../learnSPH/densities.h"
#include "../learnSPH/geometry.h"
#include "../learnSPH/io.h"
#include "../learnSPH/kernel.h"
#include "../learnSPH/sampling.h"
#include "../learnSPH/time_integration.h"

#include <chrono>

int main() {
    std::cout << "Welcome to the learnSPH framework!!" << std::endl;

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

    const double fluid_mass = 1000.0;
    const double fluid_volume = 1.0;
    const double fluid_density = fluid_mass / fluid_volume;

    // create cubic spline kernel
    learnSPH::kernel::CubicSplineKernel cubic_kernel(h);

    // Load simulation geometry
    std::vector<Eigen::Vector3d> boundary_particles_positions;
    // learnSPH::geometry::load_n_sample_boundary(boundary_particles_positions
    // 										     ,
    // "./res/box.obj"
    // , boundary_sampling_distance);

    std::cout << "Number of boundary particles" << std::endl;
    std::cout << boundary_particles_positions.size() << std::endl;

    // Sample the fluid with particles
    std::vector<Eigen::Vector3d> particles_positions;
    learnSPH::sampling::fluid_box(
            particles_positions, Eigen::Vector3d(0.0, 0.0, 0.0),
            Eigen::Vector3d(1.0, 1.0, 1.0), fluid_sampling_distance);

    int end_of_box1 = particles_positions.size() - 1;

    learnSPH::sampling::fluid_box(particles_positions, {1.1, 0.0, 0.0},
                                {2.1, 1.0, 1.0}, fluid_sampling_distance);
    int end_of_box2 = particles_positions.size() - 1;

    std::cout << "Number of fluid particles" << std::endl;
    std::cout << particles_positions.size() << std::endl;

    double fluid_particle_mass = 2 * fluid_mass / particles_positions.size();

    // Compute boundary mass
    CompactNSearch::NeighborhoodSearch nsearch(beta);

    unsigned int point_set_id_boundary =
            nsearch.add_point_set(boundary_particles_positions.front().data(),
                                                        boundary_particles_positions.size());
    unsigned int point_set_id_fluid = nsearch.add_point_set(
            particles_positions.front().data(), particles_positions.size());

    nsearch.set_active(point_set_id_boundary, point_set_id_fluid, false);
    nsearch.set_active(point_set_id_fluid, point_set_id_fluid, true);
    nsearch.set_active(point_set_id_fluid, point_set_id_boundary, true);
    nsearch.find_neighbors();

    CompactNSearch::PointSet const &ps_boundary =
            nsearch.point_set(point_set_id_boundary);
    CompactNSearch::PointSet const &ps_fluid =
            nsearch.point_set(point_set_id_fluid);

    std::vector<double> boundary_particles_masses(boundary_particles_positions.size());

    learnSPH::densities::compute_boundary_masses(
            boundary_particles_masses, boundary_particles_positions,
            point_set_id_boundary, ps_boundary, fluid_density, cubic_kernel);

    double boundary_volume = std::min(0.0,
           							 std::accumulate(boundary_particles_masses.begin(),
													boundary_particles_masses.end(),
													decltype(boundary_particles_masses)::value_type(0.0)) /
                    				fluid_density);

    double dt_cfl = 0.005;
    double dt_default = 0.001;
    double dt;
    double t_next_frame = 0;
    double t_between_frames = 0.001;
    double t_simulation = 0;
    double B = 1000 * 1.02;
    double v_f = 1.002 * 10e-6;
    double v_b = 0;
    Eigen::Vector3d gravity = Eigen::Vector3d(0.0, 0.0, 0.0);

    std::vector<double> particles_densities(particles_positions.size());
    std::vector<Eigen::Vector3d> particles_velocities(
            particles_positions.size(), Eigen::Vector3d(0.0, 0.0, 0.0));


    for (int i = 0; i <= end_of_box1; ++i) {
        particles_velocities[i] = {1, 0, 0};
    }
    for (int i = end_of_box1 + 1; i <= end_of_box2; ++i) {
        particles_velocities[i] = {-1, 0, 0};
    }


    std::vector<Eigen::Vector3d> particles_accelerations(particles_positions.size());
    std::vector<double> particles_pressure(particles_positions.size());

    learnSPH::acceleration::Acceleration acceleration(B, v_f, v_b, h, fluid_density, gravity, cubic_kernel);
    learnSPH::timeIntegration::semiImplicitEuler semImpEuler(particle_radius);

    const std::string filename = "./res/wcsph/wcsph_0.vtk";
    learnSPH::write_particles_to_vtk(filename, particles_positions);

    // Simulation loop
    while (t_simulation < 5) {

        std::cout << "t_simulation " << t_simulation << std::endl;
        std::cout << "v_max " << semImpEuler.v_max << std::endl;
        // Compute dt
        dt_cfl = 0.5 * particle_radius * (1 / semImpEuler.v_max);
        dt = std::min(dt_cfl, dt_default);
        std::cout << "dt " << dt << std::endl;

        // Find neighbors
        auto start_neighborhood = std::chrono::high_resolution_clock::now();

        nsearch.find_neighbors();

        auto stop_neighborhood = std::chrono::high_resolution_clock::now();
        auto duration_neighborhood =
                std::chrono::duration_cast<std::chrono::microseconds>(
                        stop_neighborhood - start_neighborhood);
        std::cout << "Compute nsearch time: ";
        std::cout << duration_neighborhood.count();
        std::cout << " microseconds." << std::endl;

        auto start_density = std::chrono::high_resolution_clock::now();
        // Compute fluid particles densities
        learnSPH::densities::compute_fluid_density(
                particles_densities, particles_positions, boundary_particles_positions,
                boundary_particles_masses, point_set_id_fluid, ps_fluid,
                point_set_id_boundary, ps_boundary, fluid_particle_mass, cubic_kernel);

        auto stop_density = std::chrono::high_resolution_clock::now();
        auto duration_density =
                std::chrono::duration_cast<std::chrono::microseconds>(stop_density -
                                                                                                                            start_density);
        std::cout << "Compute density time: ";
        std::cout << duration_density.count();
        std::cout << " microseconds." << std::endl;

        // Compute acceleration
        auto start_pressure = std::chrono::high_resolution_clock::now();
        acceleration.pressure(particles_pressure, particles_densities,
                                                    fluid_density);
        auto stop_pressure = std::chrono::high_resolution_clock::now();
        auto duration_pressure =
                std::chrono::duration_cast<std::chrono::microseconds>(stop_pressure -
                                                                                                                            start_pressure);
        std::cout << "Compute pressure time: ";
        std::cout << duration_pressure.count();
        std::cout << " microseconds." << std::endl;

        auto start_acceleration = std::chrono::high_resolution_clock::now();
	
        acceleration.accelerations(
                particles_accelerations, particles_densities, particles_pressure,
                point_set_id_fluid, point_set_id_boundary, ps_fluid, ps_boundary,
                particles_positions, boundary_particles_positions, particles_velocities,
                boundary_volume, fluid_particle_mass);
	
        auto stop_acceleration = std::chrono::high_resolution_clock::now();
        auto duration_acceleration =
                std::chrono::duration_cast<std::chrono::microseconds>(
                        stop_acceleration - start_acceleration);
        std::cout << "Compute acceleration time: ";
        std::cout << duration_acceleration.count();
        std::cout << " microseconds." << std::endl;

        // Integrate
        auto start_integrate = std::chrono::high_resolution_clock::now();
        semImpEuler.integrationStep(particles_positions, particles_velocities,
                                                                particles_accelerations, dt);
        auto stop_integrate = std::chrono::high_resolution_clock::now();
        auto duration_integrate =
                std::chrono::duration_cast<std::chrono::microseconds>(stop_integrate -
                                                                                                                            start_integrate);
        std::cout << "Compute integration time: ";
        std::cout << duration_integrate.count();
        std::cout << " microseconds." << std::endl;

        // Increment t
        t_simulation += dt;

        // Save output
        if (t_simulation >= t_next_frame) {
            const std::string filename = "./res/wcsph/wcsph_" + std::to_string((int)(t_simulation * 1000)) +".vtk";
            learnSPH::write_particles_to_vtk(filename, particles_positions);
            t_next_frame += t_between_frames;
        }
    }
    return 0;
}