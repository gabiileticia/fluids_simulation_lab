#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <iostream>
#include <memory>
#include <ostream>
#include <sstream>
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

int main()
{
    bool filesetup               = true;
    std::string sparse_logFile   = "Mcubes_sparse_log.txt";
    std::string dense_logFile    = "Mcubes_dense_log.txt";
    std::string setup_logFile    = "Mcubes_setup_log.txt";
    std::string levelmap_logFile = "Mcubes_levelmap_log.txt";
    std::string levelset_logFile = "Mcubes_levelset_log.txt";

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
    // sim_setup.dam_overspill();
    //  sim_setup.our_simulation_scene();
    // sim_setup.mcubes_stress_test_scene();

    double particle_diameter          = 2.0 * sim_setup.particle_radius;
    double fluid_sampling_distance    = particle_diameter;
    double boundary_sampling_distance = 0.8 * particle_diameter;
    double h                          = 1.2 * particle_diameter;
    double beta                       = 2.0 * h;
    double epsilon                    = 1e-6;
    std::vector<std::array<int, 2>> boundaries;

    using namespace learnSPH;

    // Marching cubes setup
    double c                = 0.55;
    double cell_width       = 1.25f * sim_setup.particle_radius; // ideal: 1.0 <= cell_width <= 1.25
    Eigen::Vector3d bborder = Eigen::Vector3d(1.5 * beta, 1.5 * beta, 1.5 * beta);

    Eigen::Vector3d min_fluid_reco;
    Eigen::Vector3d max_fluid_reco;

    // instantiating some classes
    std::cout << sim_setup.assignment << ", " << simulation_timestamp << "\n";
    utils::create_simulation_folder(sim_setup.assignment, simulation_timestamp);
    const std::string log_file =
        "./res/" + sim_setup.assignment + "/" + simulation_timestamp + "/log.txt";

    learnSPH::kernel::CubicSplineKernel cubic_kernel(h, beta);
    learnSPH::acceleration::Acceleration acceleration(sim_setup.B, sim_setup.v_f, sim_setup.v_b, h,
                                                      sim_setup.fluid_rest_density,
                                                      sim_setup.gravity, cubic_kernel);
    learnSPH::timeIntegration::semiImplicitEuler semImpEuler(
        sim_setup.particle_radius, sim_setup.objects,
        {sim_setup.sim_boundary_min, sim_setup.sim_boundary_max});

    // Load simulation geometry
    std::vector<Eigen::Vector3d> boundary_particles_positions;
    geometry::load_n_sample_boundary(boundary_particles_positions, sim_setup.objects, boundaries,
                                     boundary_sampling_distance);

    std::cout << "Number of boundary particles" << std::endl;
    std::cout << boundary_particles_positions.size() << std::endl;

    // Sample the fluid with particles
    std::vector<Eigen::Vector3d> particles_positions;
    std::vector<Eigen::Vector3d> particles_velocities;

    double fluid_particle_mass;

    geometry::load_n_sample_fluids(particles_positions, particles_velocities, sim_setup.fluid_begin,
                                   sim_setup.fluid_end, fluid_sampling_distance,
                                   sim_setup.fluid_velocities);

    if (sim_setup.sampleMassbyFluid) {
        double fluid_volume = (sim_setup.fluid_end[0].x() - sim_setup.fluid_begin[0].x()) *
                              (sim_setup.fluid_end[0].y() - sim_setup.fluid_begin[0].y()) *
                              (sim_setup.fluid_end[0].z() - sim_setup.fluid_begin[0].z());

        fluid_particle_mass =
            fluid_volume * sim_setup.fluid_rest_density / particles_positions.size();
    } else {
        fluid_particle_mass =
            learnSPH::utils::particle_mass(sim_setup.fluid_rest_density, fluid_sampling_distance);
    }

    std::cout << "Number of fluid particles" << std::endl;
    std::cout << particles_positions.size() << std::endl;
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

    std::vector<learnSPH::emitter::Emitter> emitters;
    std::vector<std::array<int, 3>> emit_mark;
    std::vector<bool> deleteFlag(particles_positions.size());

    if (sim_setup.emitters.size() > 0)
        for (int i = 0; i < sim_setup.emitters.size(); i++) {
            learnSPH::emitter::Emitter em(
                sim_setup.emitters[i].dir, sim_setup.emitters[i].origin, sim_setup.emitters[i].r,
                sim_setup.particle_radius, sim_setup.emitters[i].velocity,
                sim_setup.emitters[i].emit_counter, emit_mark, particles_positions,
                particles_accelerations, particles_velocities, particles_densities,
                particles_pressure, deleteFlag, point_set_id_fluid, nsearch);
            emitters.push_back(em);

            if (sim_setup.emitter_shield) {
                learnSPH::utils::create_emitter_shield(
                    em.rotation_matrix, sim_setup.emitters[i].origin, sim_setup.emitters[i].r,
                    boundary_particles_positions, sim_setup.particle_radius, point_set_id_boundary,
                    nsearch);
            }
        }

    nsearch.find_neighbors();
    // Compute boundary masses
    std::vector<double> boundary_particles_masses(boundary_particles_positions.size());
    densities::compute_boundary_masses(boundary_particles_masses, boundary_particles_positions,
                                       point_set_id_boundary, ps_boundary,
                                       sim_setup.fluid_rest_density, cubic_kernel);
    // keeping track of number of elements which will be deleted
    int count_del        = 0;
    int count_del_it_sum = 0;

    int maxSteps    = sim_setup.simTime / sim_setup.t_between_frames;
    int stepCounter = 0;

    const std::string boundary_file =
        "./res/" + sim_setup.assignment + "/" + simulation_timestamp + "/boundary_particles.vtk";

    write_particles_to_vtk(boundary_file, boundary_particles_positions);

    double hcp_z = sim_setup.particle_radius * (2 * std::sqrt(6)) / 3;

    std::ostringstream msg;

    msg << "Pressure solver: ";
    if (sim_setup.pressure_solver_method == 1) {
        msg << "position based\n";
        msg << "Iterations pbf: " << sim_setup.n_iterations_pbf << "\n";
    } else {
        msg << "weakly compressible\n";
    }

    msg << "c: " << c << "\n";
    msg << "cell width: " << cell_width << "\n";
    msg << "particle radius: " << sim_setup.particle_radius << "\n";
    // log simsetup settings
    msg << "delta t default: " << sim_setup.dt_default << "\n"
        << "frame time: " << sim_setup.t_between_frames << "\n"
        << "Rest density: " << sim_setup.B << "\n"
        << "viscosity fluid: " << sim_setup.v_f << "\n"
        << "viscosity boundaries: " << sim_setup.v_b << "\n"
        << "gravity: " << sim_setup.gravity.x() << ", " << sim_setup.gravity.y() << ", "
        << sim_setup.gravity.z() << "\n"
        << "Sim boundary active: " << (sim_setup.simbound_active ? "yes" : "no") << "\n";

    if (sim_setup.surface_tension) {
        msg << "Cohesion value: " << sim_setup.cohesion_coefficient
            << "\nAdhesion value: " << sim_setup.adhesion_coefficient << "\n";
    }

    utils::logMessage(msg.str(), log_file);

    // Simulation loop
    while (t_simulation < sim_setup.simTime) {

        // Compute dt
        dt_cfl =
            0.5 * sim_setup.particle_radius * (1 / std::min(100.0, std::sqrt(semImpEuler.v_max)));
        dt = std::min(dt_cfl, sim_setup.dt_default);

        for (int i = 0; i < emitters.size(); i++) {
            if ((t_simulation - emitters[i].last_emit) * emitters[i].emit_velocity >
                    (hcp_z * sim_setup.emitters[i].emission_freq) &&
                emitters[i].emit_counter > 0) {
                emitters[i].emit_particles_alternating(t_simulation, i);
            }
        }

        if (particles_positions.size() > 0) {
            // Find neighbors
            nsearch.find_neighbors();

            // Compute fluid particles densities
            learnSPH::densities::compute_fluid_density(
                particles_densities, particles_positions, boundary_particles_positions,
                boundary_particles_masses, point_set_id_fluid, ps_fluid, point_set_id_boundary,
                fluid_particle_mass, cubic_kernel);

            // Compute acceleration
            acceleration.pressure(particles_pressure, particles_densities,
                                  sim_setup.fluid_rest_density);

            acceleration.accelerations(
                particles_accelerations, particles_densities, particles_pressure,
                point_set_id_fluid, point_set_id_boundary, ps_fluid, ps_boundary,
                particles_positions, boundary_particles_positions, particles_velocities,
                boundary_particles_masses, sim_setup.fluid_rest_density, fluid_particle_mass);

            for (int i = 0; i < emit_mark.size(); ++i) {
                double d = (emitters[emit_mark[i][2]].dir.dot(particles_positions[emit_mark[i][0]] -
                                                              emitters[emit_mark[i][2]].origin)) /
                           emitters[emit_mark[i][2]].dir.norm();
                if (d > 3 * particle_diameter)
                    emit_mark.erase(emit_mark.begin() + i);
            }

            for (int i = 0; i < emit_mark.size(); i++) {
                for (int j = emit_mark[i][0]; j < emit_mark[i][1]; j++) {
                    particles_accelerations[j] = {0, 0, 0};
                }
            }
            // Integrate
            semImpEuler.integrationStep(particles_positions, particles_velocities,
                                        particles_accelerations, deleteFlag, dt, count_del,
                                        min_fluid_reco, max_fluid_reco);
        } else {
            dt = sim_setup.t_between_frames;
        }

        if (count_del > 0 && (sim_setup.objects.size() > 0 || sim_setup.simbound_active)) {
            learnSPH::utils::deleteOutOfBounds(particles_positions, particles_velocities,
                                               particles_accelerations, particles_densities,
                                               particles_pressure, deleteFlag, count_del);
            nsearch.resize_point_set(point_set_id_fluid, particles_positions.front().data(),
                                     particles_positions.size());
            count_del_it_sum += count_del;
        }

        // Increment t
        t_simulation += dt;

        // Save output
        if (t_simulation >= t_next_frame) {
            stepCounter++;

            const std::string filename_dense =
                "./res/" + sim_setup.assignment + "/" + simulation_timestamp + "/dense_sim_" +
                std::to_string((int)(t_simulation * 1000000)) + ".vtk";
            const std::string filename_sparse =
                "./res/" + sim_setup.assignment + "/" + simulation_timestamp + "/sparse_sim_" +
                std::to_string((int)(t_simulation * CLOCKS_PER_SEC)) + ".vtk";
            const std::string fileprefix =
                "./res/" + sim_setup.assignment + "/" + simulation_timestamp + "/";
            const std::string particles_filename =
                "./res/" + sim_setup.assignment + "/" + simulation_timestamp + "/particles_" +
                std::to_string((int)(t_simulation * 1000000)) + ".vtk";

            uint nx = ((max_fluid_reco.x() + bborder.x()) - (min_fluid_reco.x() - bborder.x())) /
                          cell_width +
                      1;
            uint ny = ((max_fluid_reco.y() + bborder.y()) - (min_fluid_reco.y() - bborder.y())) /
                          cell_width +
                      1;
            uint nz = ((max_fluid_reco.z() + bborder.z()) - (min_fluid_reco.z() - bborder.z())) /
                          cell_width +
                      1;

            nsearch.find_neighbors();

            learnSPH::densities::compute_fluid_density_surface_reco(
                fluid_densities_for_surface_reco, particles_positions, point_set_id_fluid, ps_fluid,
                cubic_kernel);

            // sparse_logFile = fileprefix + sparse_logFile;

            // setup for mcubes - two separate objects so no sideeffects
            // need to get rid of logging here eventually
            auto start_setup = std::chrono::high_resolution_clock::now();
            std::unordered_map<uint64_t, double> level_map;
            std::vector<double> level_set((nx + 1) * (ny + 1) * (nz + 1), -c);
            learnSPH::theta_functions::FluidThetaFunction fluidSDF(cubic_kernel, c, cell_width,
                                                                   beta, nx + 1, ny + 1, nz + 1);
            learnSPH::surface::MarchingCubes mcubes_sparse(cell_width, nx, ny, nz,
                                                           min_fluid_reco - bborder, epsilon);

            learnSPH::surface::MarchingCubes mcubes_dense(cell_width, nx, ny, nz,
                                                          min_fluid_reco - bborder, epsilon);
            auto end_setup = std::chrono::high_resolution_clock::now();

            // // creating level set vector
            auto start_levelset = std::chrono::high_resolution_clock::now();
            fluidSDF.computeLevelSet(level_set, particles_positions,
                                     fluid_densities_for_surface_reco, min_fluid_reco - bborder);
            auto end_levelset = std::chrono::high_resolution_clock::now();

            // creating level set hash map
            auto start_levelmap = std::chrono::high_resolution_clock::now();
            fluidSDF.computeLevelMap(level_map, particles_positions,
                                     fluid_densities_for_surface_reco, min_fluid_reco - bborder);
            auto end_levelmap = std::chrono::high_resolution_clock::now();

            // computing isosurface dense
            auto start_dense = std::chrono::high_resolution_clock::now();
            mcubes_dense.get_Isosurface(level_set);
            std::ostringstream dense_infix;
            dense_infix << mcubes_dense.intersections.size() << ","
                        << mcubes_dense.triangles.size();
            auto end_dense = std::chrono::high_resolution_clock::now();

            // computing isosurface sparse
            auto start_sparse = std::chrono::high_resolution_clock::now();
            mcubes_sparse.get_Isosurface_sparse(level_map);
            std::ostringstream sparse_infix;
            sparse_infix << mcubes_sparse.intersections.size() << ","
                         << mcubes_sparse.triangles.size();
            auto end_sparse = std::chrono::high_resolution_clock::now();

            // computing normals - logging maybe interesting but i don't expect differences here
            mcubes_sparse.compute_normals();
            mcubes_dense.compute_normals();

            std::chrono::duration<double> duration_setup    = end_setup - start_setup;
            std::chrono::duration<double> duration_dense    = end_dense - start_dense;
            std::chrono::duration<double> duration_sparse   = end_sparse - start_sparse;
            std::chrono::duration<double> duration_levelmap = end_levelmap - start_levelmap;
            std::chrono::duration<double> duration_leveset  = end_levelset - start_levelset;

            std::ostringstream msg_setup;
            std::ostringstream msg_dense;
            std::ostringstream msg_sparse;
            std::ostringstream msg_levelmap;
            std::ostringstream msg_levelset;

            using namespace learnSPH::utils;

            if (filesetup) {
                std::string header =
                    "#Particles,#Grid Vertices,#Mesh Vertices,#Mesh Triangles,Duration in s";
                std::string sparse_header = "#Level Map," + header;
                logMessage(sparse_header, fileprefix + sparse_logFile);
                logMessage(header, fileprefix + dense_logFile);
                logMessage(header, fileprefix + setup_logFile);
                logMessage(header, fileprefix + levelmap_logFile);
                logMessage(header, fileprefix + levelset_logFile);
                filesetup = false;
            }

            // information same for every entry
            // save levelmap size for sparse data for more compareability to dense -> large map
            //      leads to sparse being slower then dense version!
            std::ostringstream prefix;
            prefix << particles_positions.size() << "," << level_set.size();
            msg_sparse << level_map.size() << ',' << prefix.str() << "," << sparse_infix.str()
                       << "," << duration_sparse.count();
            msg_dense << prefix.str() << "," << dense_infix.str() << "," << duration_dense.count();
            msg_setup << prefix.str() << ","
                      << ","
                      << "," << duration_setup.count();
            msg_levelmap << prefix.str() << ","
                         << ","
                         << "," << duration_levelmap.count();
            msg_levelset << prefix.str() << ","
                         << ","
                         << "," << duration_leveset.count();

            logMessage(msg_dense.str(), fileprefix + dense_logFile);
            logMessage(msg_sparse.str(), fileprefix + sparse_logFile);
            logMessage(msg_setup.str(), fileprefix + setup_logFile);
            logMessage(msg_levelmap.str(), fileprefix + levelmap_logFile);
            logMessage(msg_levelset.str(), fileprefix + levelset_logFile);

            //
            // write_tri_mesh_to_vtk(filename_dense, mcubes_dense.intersections,
            //                       mcubes_dense.triangles, mcubes_dense.intersectionNormals);
            // write_tri_mesh_to_vtk(filename_sparse, mcubes_sparse.intersections,
            //                       mcubes_sparse.triangles, mcubes_sparse.intersectionNormals);
            // write_particles_to_vtk(particles_filename, particles_positions, particles_densities,
            //                        particles_velocities);
            t_next_frame += sim_setup.t_between_frames;

            if (count_del_it_sum > 0) {
                std::cout << "Deleting " << count_del_it_sum << " elements from particle vectors."
                          << std::endl;
                count_del_it_sum = 0;
            }
            auto progress_msg = utils::updateProgressBar(stepCounter, maxSteps, 75);
            std::cout << progress_msg.str() << std::endl;
            logMessage(progress_msg.str(), log_file);
        }
    }
    return 0;
}