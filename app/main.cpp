#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <memory>
#include <ostream>
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include "../extern/CompactNSearch/include/CompactNSearch/CompactNSearch"
#include "../learnSPH/acceleration.h"
#include "../learnSPH/densities.h"
#include "../learnSPH/geometry.h"
#include "../learnSPH/io.h"
#include "../learnSPH/kernel.h"
#include "../learnSPH/marching_cubes.h"
#include "../learnSPH/sampling.h"
#include "../learnSPH/simulations_setup.h"
#include "../learnSPH/theta_functions.h"
#include "../learnSPH/time_integration.h"
#include "../learnSPH/utils.h"
#include "../learnSPH/theta_functions.h"
#include "../learnSPH/marching_cubes.h"

int main()
{

    std::cout << "Welcome to the learnSPH framework!!" << std::endl;

    double dt_cfl, dt;
    double t_next_frame = 0;
    double t_simulation = 0;
    std::string simulation_timestamp;

    // Setting up simulation
    learnSPH::simulations_setup::Simulations sim_setup;
    //     sim_setup.simple_cube();
    //     sim_setup.simple_cube_with_fluid_viscosity();
    //     sim_setup.cubes_colision();
    //     sim_setup.just_gravity();
    //     sim_setup.gravity_with_floor();
    //     sim_setup.gravity_with_floor_boundary_viscosity();
        sim_setup.dam_break();
    // sim_setup.our_simulation_scene();

    double particle_diameter          = 2.0 * sim_setup.particle_radius;
    double fluid_sampling_distance    = particle_diameter;
    double boundary_sampling_distance = 0.8 * particle_diameter;
    double h                          = 1.2 * particle_diameter;
    double beta                       = 2.0 * h;

    double fluid_volume = 0.0;
    for (int i = 0; i < sim_setup.fluid_begin.size(); ++i) {
        fluid_volume += (sim_setup.fluid_end[i].x() - sim_setup.fluid_begin[i].x()) *
                        (sim_setup.fluid_end[i].y() - sim_setup.fluid_begin[i].y()) *
                        (sim_setup.fluid_end[i].z() - sim_setup.fluid_begin[i].z());
    }
    double fluid_mass = fluid_volume * sim_setup.fluid_rest_density;

    using namespace learnSPH;

    // Marching cubes setup
    double c = 0.55;
    double cell_width = 1.25 * sim_setup.particle_radius;
    Eigen::Vector3d boundary_border = Eigen::Vector3d(0.2, 0.2, 0.2);
    
    uint nx = ((sim_setup.boundaries[0].max.x() + boundary_border.x()) - (sim_setup.boundaries[0].min.x() - boundary_border.x())) / cell_width + 1;
    uint ny = ((sim_setup.boundaries[0].max.y() + boundary_border.y()) - (sim_setup.boundaries[0].min.y() - boundary_border.y())) / cell_width + 1;
    uint nz = ((sim_setup.boundaries[0].max.z() + boundary_border.z()) - (sim_setup.boundaries[0].min.z() - boundary_border.z())) / cell_width + 1;

    // instantiating some classes
    std::cout << sim_setup.assignment << ", " << simulation_timestamp << "\n";
    utils::create_simulation_folder(sim_setup.assignment, simulation_timestamp);

    learnSPH::kernel::CubicSplineKernel cubic_kernel(h);
    learnSPH::acceleration::Acceleration acceleration(sim_setup.B, sim_setup.v_f, sim_setup.v_b, h,
                                                      sim_setup.fluid_rest_density,
                                                      sim_setup.gravity, cubic_kernel);
    learnSPH::timeIntegration::semiImplicitEuler semImpEuler(sim_setup.particle_radius,
                                                             sim_setup.boundaries);

    

    // Load simulation geometry
    std::vector<Eigen::Vector3d> boundary_particles_positions;
    geometry::load_n_sample_boundary(boundary_particles_positions, sim_setup.boundaries,
                                     boundary_sampling_distance);

    std::cout << "Number of boundary particles" << std::endl;
    std::cout << boundary_particles_positions.size() << std::endl;

    // Sample the fluid with particles
    std::vector<Eigen::Vector3d> particles_positions;
    std::vector<Eigen::Vector3d> particles_velocities;

    geometry::load_n_sample_fluids(particles_positions, particles_velocities, sim_setup.fluid_begin,
                                   sim_setup.fluid_end, fluid_sampling_distance,
                                   sim_setup.fluid_velocities);

    std::cout << "Number of fluid particles" << std::endl;
    std::cout << particles_positions.size() << std::endl;

    double fluid_particle_mass = fluid_mass / particles_positions.size();
    std::cout << "fluid_particle_mass: " << fluid_particle_mass << std::endl;

    std::vector<double> particles_densities(particles_positions.size());
    std::vector<double> fluid_densities_for_surface_reco(particles_positions.size());
    std::vector<Eigen::Vector3d> particles_accelerations(particles_positions.size());
    std::vector<double> particles_pressure(particles_positions.size());

    // Setting up neighborhood search
    CompactNSearch::NeighborhoodSearch nsearch(beta);
    unsigned int point_set_id_boundary = nsearch.add_point_set(
        boundary_particles_positions.front().data(), boundary_particles_positions.size());

    unsigned int point_set_id_fluid =
        nsearch.add_point_set(particles_positions.front().data(), particles_positions.size());

    nsearch.set_active(point_set_id_boundary, point_set_id_fluid, false);
    nsearch.set_active(point_set_id_boundary, point_set_id_boundary, true);
    nsearch.set_active(point_set_id_fluid, point_set_id_fluid, true);
    nsearch.set_active(point_set_id_fluid, point_set_id_boundary, true);
    nsearch.find_neighbors();

    CompactNSearch::PointSet const &ps_boundary = nsearch.point_set(point_set_id_boundary);
    CompactNSearch::PointSet const &ps_fluid    = nsearch.point_set(point_set_id_fluid);

    // Compute boundary masses
    std::vector<double> boundary_particles_masses(boundary_particles_positions.size());
    densities::compute_boundary_masses(boundary_particles_masses, boundary_particles_positions,
                                       point_set_id_boundary, ps_boundary,
                                       sim_setup.fluid_rest_density, cubic_kernel);

    // setting up surface reconstruction
    double c                = 0.55;
    double cell_width       = 1.25 * sim_setup.particle_radius;
    Eigen::Vector3d bborder = Eigen::Vector3d(5, 5, 5) * cell_width;

    uint nx = ((sim_setup.boundaries[0].max.x() + bborder.x()) -
               (sim_setup.boundaries[0].min.x() - bborder.x())) /
                  cell_width +
              1;
    uint ny = ((sim_setup.boundaries[0].max.y() + bborder.y()) -
               (sim_setup.boundaries[0].min.y() - bborder.y())) /
                  cell_width +
              1;
    uint nz = ((sim_setup.boundaries[0].max.z() + bborder.z()) -
               (sim_setup.boundaries[0].min.z() - bborder.z())) /
                  cell_width +
              1;

    std::vector<double> level_set((nx + 1) * (ny + 1) * (nz + 1), -c);

    // keeping track of number of elements which will be deleted
    int count_del = 0;
    std::vector<bool> deleteFlag(particles_positions.size());

    int maxSteps    = 5 / sim_setup.t_between_frames;
    int stepCounter = 0;

    // Simulation loop
    while (t_simulation < 5) {

        // Compute dt
        dt_cfl =
            0.5 * sim_setup.particle_radius * (1 / std::min(100.0, std::sqrt(semImpEuler.v_max)));
        dt = std::min(dt_cfl, sim_setup.dt_default);

        // Find neighbors
        nsearch.find_neighbors();

        // Compute fluid particles densities
        learnSPH::densities::compute_fluid_density(fluid_densities_for_surface_reco,
            particles_densities, particles_positions, boundary_particles_positions,
            boundary_particles_masses, point_set_id_fluid, ps_fluid, point_set_id_boundary,
            ps_boundary, fluid_particle_mass, cubic_kernel);

        std::vector<double> level_set((nx + 1) * (ny + 1) * (nz + 1), -c);
        learnSPH::theta_functions::FluidThetaFunction fluid_theta(c, cell_width, beta, nx + 1, ny + 1, nz + 1);
        learnSPH::surface::MarchingCubes mcubes(cell_width, nx, ny, nz, sim_setup.boundaries[0].min - boundary_border, 1e-6, true);
    
        fluid_theta.compute_fluid_reconstruction(level_set, particles_positions, fluid_densities_for_surface_reco,
                sim_setup.boundaries[0].min - boundary_border, sim_setup.boundaries[0].max + boundary_border, cubic_kernel);
        mcubes.get_isosurface(level_set);
        mcubes.compute_normals(level_set);

        // Compute acceleration
        acceleration.pressure(particles_pressure, particles_densities,
                              sim_setup.fluid_rest_density);

        acceleration.accelerations(particles_accelerations, particles_densities, particles_pressure,
                                   point_set_id_fluid, point_set_id_boundary, ps_fluid, ps_boundary,
                                   particles_positions, boundary_particles_positions,
                                   particles_velocities, boundary_particles_masses,
                                   sim_setup.fluid_rest_density, fluid_particle_mass);

        // Integrate
        semImpEuler.integrationStep(particles_positions, particles_velocities,
                                    particles_accelerations, deleteFlag, dt, count_del);

        if (count_del > 0 && sim_setup.boundaries.size() > 0) {
            learnSPH::utils::deleteOutOfBounds(particles_positions, particles_velocities,
                                               particles_accelerations, particles_densities,
                                               fluid_densities_for_surface_reco,
                                               particles_pressure, deleteFlag, count_del);
            nsearch.resize_point_set(point_set_id_fluid, particles_positions.front().data(),
                                     particles_positions.size());
        }

        // Increment t
        t_simulation += dt;

        // Save output
        if (t_simulation >= t_next_frame) {
            stepCounter++;

            const std::string filename = "./res/" + sim_setup.assignment + "/" +
                                         simulation_timestamp + "/sim_" +
                                         std::to_string((int)(t_simulation * 1000000)) + ".vtk";
            // std::cout << "num pos: " << particles_positions.size() << "\n";
            // std::cout << "t_simulation: " << t_simulation << "\n t_nextframe: " << t_next_frame
            // << "\n";

            fluidSDF.computeLevelMap(level_map, particles_positions, particles_densities,
                                     sim_setup.boundaries[0].min - bborder);
            // fluidSDF.computeLevelSet(level_set, particles_positions, particles_densities,
            //   sim_setup.boundaries[0].min - bborder);
            // mcubes.get_Isosurface(level_set);
            mcubes.get_Isosurface_sparse(level_map);

            mcubes.compute_normals(particles_positions, particles_densities,
                                   sim_setup.boundaries[0].min - bborder);

            write_tri_mesh_to_vtk(filename, mcubes.intersections, mcubes.triangles,
                                  mcubes.intersectionNormals);

            //write_particles_to_vtk(filename, mcubes.debug);

            // write_particles_to_vtk(filename, particles_positions, particles_densities,
            //                        particles_velocities);
            const std::string filename = "./res/surface2/sim_" + std::to_string((int)(t_simulation * 1000000)) + ".vtk";

            learnSPH::write_tri_mesh_to_vtk(filename, mcubes.intersections, mcubes.triangles, mcubes.intersectionNormals);

            // const std::string filename = "./res/" + sim_setup.assignment + "/" +
            //                              simulation_timestamp + "/sim_" +
            //                              std::to_string((int)(t_simulation * 1000000)) + ".vtk";

            // learnSPH::write_particles_to_vtk(filename, particles_positions, particles_densities,
            //                                  particles_velocities);
            t_next_frame += sim_setup.t_between_frames;

            // level_set.clear();
            // mcubes.intersections.clear();
            // mcubes.triangles.clear();
            // mcubes.intersectionNormals.clear();

            utils::updateProgressBar(stepCounter, maxSteps, 75);
        }
    }
    return 0;
}