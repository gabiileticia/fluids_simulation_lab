#include <Eigen/Dense>
#include <chrono>
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>

#include "../extern/CompactNSearch/include/CompactNSearch/CompactNSearch"
#include "../learnSPH/acceleration.h"
#include "../learnSPH/densities.h"
#include "../learnSPH/geometry.h"
#include "../learnSPH/io.h"
#include "../learnSPH/kernel.h"
#include "../learnSPH/sampling.h"
#include "../learnSPH/time_integration.h"
#include "../learnSPH/utils.h"

int main()
{
    // declaring variables beforehand to make assignments less cluttered
    unsigned int point_set_id_boundary, point_set_id_fluid, size_particles, size_boundary;
    int count_del;
    double particle_radius, a, b, c, dt_cfl, dt, dt_default, t_next_frame, t_between_frames,
        t_simulation, B, v_f, v_b, fluid_particle_mass;
    Eigen::Vector3d fluid_begin, fluid_end, boundary_begin, boundary_end, gravity;
    std::vector<Eigen::Vector3d> particles_positions, particles_accelerations, particles_velocities,
        boundary_particle_positions;
    std::vector<double> particles_densities, particles_pressure, boundary_particles_densities,
        boundary_particles_masses;
    std::vector<bool> deleteFlags;
    std::string filename;

    // particles settings
    particle_radius                         = 0.005;
    const double particle_diameter          = 2 * particle_radius;
    const double fluid_sampling_distance    = particle_diameter;
    const double boundary_sampling_distance = 0.8 * particle_diameter;
    // smoothing length
    const double h = 1.2 * particle_diameter;
    // compact support
    const double beta = 2.0 * h;

    // Fluid settings
    const double fluid_rest_density = 1000.0;
    fluid_begin                     = {0.05, 0.05, 0.05};
    fluid_end                       = {.1, .25, .5};
    const double fluid_volume       = (fluid_end.x() - fluid_begin.x()) *
                                (fluid_end.y() - fluid_begin.y()) *
                                (fluid_end.z() - fluid_begin.z());
    const double fluid_mass = fluid_volume * fluid_rest_density;

    const std::string boundary_file = "./res/boundary.obj";
    a                               = .15;
    b                               = .8;
    c                               = 1.0;
    boundary_begin                  = {.0, .0, .0};
    boundary_end                    = {a, b, c};
    const double particle_volume =
        4.0 / 3.0 * learnSPH::kernel::PI * particle_radius * particle_radius * particle_radius;
    const double boundary_volume = 2 * particle_volume * (a * b + a * c + b * c);
    const double boundary_mass   = boundary_volume * fluid_rest_density;

    // simulation parameter
    dt_default       = 0.00025;
    t_next_frame     = 0;
    t_between_frames = 0.0005;
    t_simulation     = 0;
    B                = 1000 * 1.02;
    v_f              = 0.0025;
    v_b              = 0.0;
    gravity          = {0, 0, -9.81};

    // instantiating classes for kernel, acceleration and time integration
    using namespace learnSPH;
    kernel::CubicSplineKernel cubic_kernel(h);
    acceleration::Acceleration acceleration(B, v_f, v_b, h, fluid_rest_density, gravity,
                                            cubic_kernel);
    timeIntegration::semiImplicitEuler semImpEuler(particle_radius, boundary_end, boundary_begin);

    // Load simulation geometry
    geometry::load_n_sample_boundary(boundary_particle_positions, boundary_file,
                                     boundary_sampling_distance);
    sampling::fluid_box(particles_positions, fluid_begin, fluid_end, fluid_sampling_distance);

    size_particles = particles_positions.size();
    size_boundary  = boundary_particle_positions.size();

    // assigning correct sizes to fluid related vectors
    fluid_particle_mass = fluid_mass / size_particles;
    particles_densities.resize(size_particles);
    particles_pressure.resize(size_particles);
    particles_accelerations.resize(size_particles);
    particles_velocities.resize(size_particles, {0.0, 0.0, 0.0});
    deleteFlags.resize(size_particles, false);

    // setting up neighborhood search
    CompactNSearch::NeighborhoodSearch nsearch(beta);
    point_set_id_boundary =
        nsearch.add_point_set(boundary_particle_positions.front().data(), size_boundary);
    point_set_id_fluid = nsearch.add_point_set(particles_positions.front().data(), size_particles);

    nsearch.set_active(point_set_id_boundary, point_set_id_fluid, false);
    nsearch.set_active(point_set_id_boundary, point_set_id_boundary, true);
    nsearch.set_active(point_set_id_fluid, point_set_id_fluid, true);
    nsearch.set_active(point_set_id_fluid, point_set_id_boundary, true);
    nsearch.find_neighbors();

    CompactNSearch::PointSet const &ps_boundary = nsearch.point_set(point_set_id_boundary);
    CompactNSearch::PointSet const &ps_fluid    = nsearch.point_set(point_set_id_fluid);

    // assigning boundary masses
    const double boundary_particles_mass   = boundary_mass / size_boundary;
    const double boundary_particle_density = fluid_rest_density / size_boundary;

    // assigning correct sizes to boundary related vectors
    boundary_particles_densities.resize(boundary_particle_positions.size());
    boundary_particles_masses.resize(boundary_particle_positions.size());

    densities::compute_boundary_masses(boundary_particles_masses, boundary_particle_positions,
                                       point_set_id_boundary, ps_boundary, fluid_rest_density,
                                       cubic_kernel);

    {
        std::cout << "Welcome to the learnSPH framework!!" << std::endl;
        std::cout << "Yor current setup is:" << std::endl;
        std::cout << "particle_radius: " << particle_radius << std::endl;
        std::cout << "particle_diameter: " << particle_diameter << std::endl;
        std::cout << "fluid_sampling_distance: " << fluid_sampling_distance << std::endl;
        std::cout << "boundary_sampling_distance: " << boundary_sampling_distance << std::endl;
        std::cout << "fluid_density: " << fluid_rest_density << std::endl;
        std::cout << "fluid_begin: " << fluid_begin.transpose() << std::endl;
        std::cout << "fluid_end: " << fluid_end.transpose() << std::endl;
        std::cout << "fluid_volume: " << fluid_volume << std::endl;
        std::cout << "fluid_mass: " << fluid_mass << std::endl;
        std::cout << "boundary_file: " << boundary_file << std::endl;
        std::cout << "boundary_begin: " << boundary_begin.transpose() << std::endl;
        std::cout << "boundary_end: " << boundary_end.transpose() << std::endl;
        std::cout << "boundary_volume: " << boundary_volume << std::endl;
        std::cout << "boundary_mass: " << boundary_mass << std::endl;
        std::cout << "dt_default: " << dt_default << std::endl;
        std::cout << "B: " << B << std::endl;
        std::cout << "v_f: " << v_f << std::endl;
        std::cout << "v_b: " << v_b << std::endl;
        std::cout << "gravity: " << gravity.transpose() << std::endl;
        std::cout << "Number of boundary particles" << std::endl;
        std::cout << size_boundary << std::endl;
        std::cout << "Number of fluid particles" << std::endl;
        std::cout << size_particles << std::endl;
        std::cout << "fluid_particle_mass: " << fluid_particle_mass << std::endl;
        std::cout << "boundary_particles_mass: " << boundary_particles_mass << std::endl;
    }

    while (t_simulation < 5) {
        std::cout << "t_simulation " << t_simulation << std::endl;
        std::cout << "v_max " << semImpEuler.v_max << std::endl;
        // compute dt
        dt_cfl = 0.5 * particle_radius * (1 / std::min(50.0, semImpEuler.v_max));
        dt     = std::min(dt_cfl, dt_default);

        // Find neighbors
        nsearch.find_neighbors();

        // Compute fluid particles densities

        densities::compute_fluid_density(particles_densities, particles_positions,
                                         boundary_particle_positions, boundary_particles_masses,
                                         point_set_id_fluid, ps_fluid, point_set_id_boundary,
                                         ps_boundary, fluid_particle_mass, cubic_kernel);

        // compute acceleration
        acceleration.pressure(particles_pressure, particles_densities, fluid_rest_density);

        acceleration.accelerations(particles_accelerations, particles_densities, particles_pressure,
                                   point_set_id_fluid, point_set_id_boundary, ps_fluid, ps_boundary,
                                   particles_positions, boundary_particle_positions,
                                   particles_velocities, boundary_particles_masses,
                                   boundary_particle_density, fluid_particle_mass);

        // Integrate time step
        semImpEuler.integrationStep(particles_positions, particles_velocities,
                                    particles_accelerations, deleteFlags, dt, count_del);
        utils::deleteOutOfBounds(particles_positions, particles_velocities, particles_accelerations,
                                 particles_densities, particles_pressure, deleteFlags, count_del);
        nsearch.resize_point_set(point_set_id_fluid, particles_positions.front().data(),
                                 particles_positions.size());

        auto stop_integrate = std::chrono::high_resolution_clock::now();

        // Increment t
        t_simulation += dt;

        // save output

        if (t_simulation >= t_next_frame) {
            filename = "./res/wcsph_dam_break/wcsph_" +
                       std::to_string((int)(t_simulation * CLOCKS_PER_SEC)) + ".vtk";

            write_particles_to_vtk(filename, particles_positions, particles_pressure,
                                   particles_velocities);
            t_next_frame += t_between_frames;
        }
    }

    return 0;
}