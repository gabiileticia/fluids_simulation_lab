#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
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
#include "../learnSPH/pbf.h"
#include "../learnSPH/sampling.h"
#include "../learnSPH/simulations_setup.h"
#include "../learnSPH/theta_functions.h"
#include "../learnSPH/time_integration.h"
#include "../learnSPH/utils.h"

int main(int argc, char **argv)
{
    
    if (argc != 2) {
        std::cout << "ERROR: The program accepts only one argument! received " << argc
                  << " instead!";
        exit(-1);
    }
    
    int foo_select = std::stoi(argv[1]);

    learnSPH::simulations_setup::Simulations sim_setup;

    switch (foo_select) {
    case 0:
        sim_setup.simple_cube();
        break;
    case 1:
        sim_setup.simple_cube_with_fluid_viscosity();
        break;
    case 2:
        sim_setup.cubes_colision();
        break;
    case 3:
        sim_setup.just_gravity();
        break;
    case 4:
        sim_setup.gravity_with_floor();
        break;
    case 5:
        sim_setup.gravity_with_floor_boundary_viscosity();
        break;
    case 6:
        sim_setup.dam_break();
        break;
    case 7:
        sim_setup.our_simulation_scene();
        break;
    case 8:
        sim_setup.slope_ramp_wall_vessel();
        break;
    default:
        std::cout << "Selected undefined function index. Closing program.";
        exit(-1);
    }

    std::cout << "Welcome to the learnSPH framework!!" << std::endl;

    double dt_cfl, dt;
    double t_next_frame = 0;
    double t_simulation = 0;
    std::string simulation_timestamp;

    double particle_diameter          = 2.0 * sim_setup.particle_radius;
    double fluid_sampling_distance    = particle_diameter;
    double boundary_sampling_distance = 0.8 * particle_diameter;
    double h                          = 1.2 * particle_diameter;
    double beta                       = 2.0 * h;
    double epsilon                    = 1e-6;

    double fluid_volume = 0.0;
    for (int i = 0; i < sim_setup.fluid_begin.size(); ++i) {
        fluid_volume += (sim_setup.fluid_end[i].x() - sim_setup.fluid_begin[i].x()) *
                        (sim_setup.fluid_end[i].y() - sim_setup.fluid_begin[i].y()) *
                        (sim_setup.fluid_end[i].z() - sim_setup.fluid_begin[i].z());
    }
    double fluid_mass = fluid_volume * sim_setup.fluid_rest_density;

    using namespace learnSPH;

    std::ostringstream msg;

    // Marching cubes setup
    int surface_reco_method = 1; // 0:dense 1:sparse
    double c                = 0.55;
    double cell_width       = 1.25 * sim_setup.particle_radius;
    Eigen::Vector3d bborder = Eigen::Vector3d(1.5 * beta, 1.5 * beta, 1.5 * beta);

    msg << "Mcubes version: ";
    msg << ((surface_reco_method == 0) ? "dense" : "sparse") << "\n";

    Eigen::Vector3d min_fluid_reco;
    Eigen::Vector3d max_fluid_reco;

    std::vector<double> fluid_densities_for_surface_reco;

    // PBF setup
    int pressure_solver_method = 1; // 0:wcsph 1:pbf
    int n_iteractions_pbf      = 5;
    std::vector<double> pbf_s;
    std::vector<double> pbf_c;
    std::vector<double> pbf_lambda;
    std::vector<Eigen::Vector3d> pbf_dx;
    std::vector<Eigen::Vector3d> last_particles_positions;

    msg << "Pressure solver: ";
    msg << ((pressure_solver_method == 0) ? "weakly compressible" : "positions based") << "\n";
    if (pressure_solver_method == 1) {
        msg << "Iterations pbf: " << n_iteractions_pbf << "\n";
    }

    msg << "c: " << c << "\n";
    msg << "cell width: " << cell_width << "\n";
    // log simsetup settings
    msg << "delta t default: " << sim_setup.dt_default << "\n";
    msg << "frame time: " << sim_setup.t_between_frames << "\n";
    msg << "Rest density: " << sim_setup.B << "\n",
        msg << "viscosity fluid: " << sim_setup.v_f << "\n";
    msg << "viscosity boundaries: " << sim_setup.v_b << "\n";
    msg << "gravity: " << sim_setup.gravity.x() << ", " << sim_setup.gravity.y() << ", "
        << sim_setup.gravity.z() << "\n";
    msg << "Sim boundary active: " << (sim_setup.simbound_active ? "yes" : "no") << "\n";

    // instantiating some classes
    utils::create_simulation_folder(sim_setup.assignment, simulation_timestamp);

    const std::string log_file =
        "./res/" + sim_setup.assignment + "/" + simulation_timestamp + "/log.txt";

    kernel::CubicSplineKernel cubic_kernel(h);
    acceleration::Acceleration acceleration(sim_setup.B, sim_setup.v_f, sim_setup.v_b, h,
                                            sim_setup.fluid_rest_density, sim_setup.gravity,
                                            cubic_kernel);
    timeIntegration::semiImplicitEuler semImpEuler(
        sim_setup.particle_radius, sim_setup.objects,
        {sim_setup.sim_boundary_min, sim_setup.sim_boundary_max});

    // Load simulation geometry
    std::vector<Eigen::Vector3d> boundary_particles_positions;
    geometry::load_n_sample_boundary(boundary_particles_positions, sim_setup.objects,
                                     boundary_sampling_distance);

    msg << "Number of boundary particles"
        << "\n";
    msg << boundary_particles_positions.size() << "\n";

    // Sample the fluid with particles
    std::vector<Eigen::Vector3d> particles_positions;
    std::vector<Eigen::Vector3d> particles_velocities;

    geometry::load_n_sample_fluids(particles_positions, particles_velocities, sim_setup.fluid_begin,
                                   sim_setup.fluid_end, fluid_sampling_distance,
                                   sim_setup.fluid_velocities);

    msg << "Number of fluid particles"
        << "\n";
    msg << particles_positions.size() << "\n";

    double fluid_particle_mass = fluid_mass / particles_positions.size();
    msg << "fluid particles mass: " << fluid_particle_mass << "\n";

    std::vector<double> particles_densities(particles_positions.size());
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
    // keeping track of number of elements which will be deleted
    int count_del = 0;
    std::vector<bool> deleteFlag(particles_positions.size());

    int maxSteps    = 5 / sim_setup.t_between_frames;
    int stepCounter = 0;

    utils::logMessage(msg.str(), log_file);
    std::cout << msg.str();

    nsearch.find_neighbors();
    // Simulation loop
    while (t_simulation < 5) {

        // Compute dt
        dt_cfl =
            0.5 * sim_setup.particle_radius * (1 / std::min(100.0, std::sqrt(semImpEuler.v_max)));
        dt = std::min(dt_cfl, sim_setup.dt_default);

        // Compute fluid particles densities
        learnSPH::densities::compute_fluid_density(
            particles_densities, particles_positions, boundary_particles_positions,
            boundary_particles_masses, point_set_id_fluid, ps_fluid, point_set_id_boundary,
            fluid_particle_mass, cubic_kernel);

        // Compute acceleration
        if (pressure_solver_method == 0) {
            acceleration.pressure(particles_pressure, particles_densities,
                                  sim_setup.fluid_rest_density);

            acceleration.accelerations(
                particles_accelerations, particles_densities, particles_pressure,
                point_set_id_fluid, point_set_id_boundary, ps_fluid, ps_boundary,
                particles_positions, boundary_particles_positions, particles_velocities,
                boundary_particles_masses, sim_setup.fluid_rest_density, fluid_particle_mass);
        } else if (pressure_solver_method == 1) {
            acceleration.pbf_accelerations(particles_accelerations, particles_densities,
                                           point_set_id_fluid, point_set_id_boundary, ps_fluid,
                                           particles_positions, boundary_particles_positions,
                                           particles_velocities, boundary_particles_masses,
                                           sim_setup.fluid_rest_density, fluid_particle_mass);
            last_particles_positions = particles_positions;
        }

        // Integrate
        semImpEuler.integrationStep(particles_positions, particles_velocities,
                                    particles_accelerations, deleteFlag, dt, count_del,
                                    min_fluid_reco, max_fluid_reco);

        // Find neighbors
        nsearch.find_neighbors();

        if (pressure_solver_method == 1) {
            for (int i = 0; i < n_iteractions_pbf; i++) {
                learnSPH::densities::compute_fluid_density(
                    particles_densities, particles_positions, boundary_particles_positions,
                    boundary_particles_masses, point_set_id_fluid, ps_fluid, point_set_id_boundary,
                    fluid_particle_mass, cubic_kernel);

                learnSPH::pbf::compute_c(pbf_c, particles_densities, sim_setup.fluid_rest_density);

                learnSPH::pbf::compute_s(pbf_s, particles_positions, boundary_particles_positions,
                                         fluid_particle_mass, sim_setup.fluid_rest_density,
                                         boundary_particles_masses, sim_setup.fluid_rest_density,
                                         cubic_kernel, point_set_id_fluid, ps_fluid,
                                         point_set_id_boundary);

                learnSPH::pbf::compute_lambda(pbf_lambda, pbf_c, pbf_s, epsilon);

                learnSPH::pbf::compute_dx(pbf_dx, sim_setup.fluid_rest_density, fluid_particle_mass,
                                          pbf_lambda, cubic_kernel, boundary_particles_masses,
                                          point_set_id_fluid, ps_fluid, point_set_id_boundary,
                                          particles_positions, boundary_particles_positions);

                learnSPH::pbf::update_positions(particles_positions, pbf_dx);
            }

            learnSPH::pbf::update_velocities(particles_positions, last_particles_positions,
                                             particles_velocities, dt);
        }
        if (count_del > 0 && (sim_setup.objects.size() > 0 || sim_setup.simbound_active)) {
            learnSPH::utils::deleteOutOfBounds(particles_positions, particles_velocities,
                                               particles_accelerations, particles_densities,
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

            uint nx = ((max_fluid_reco.x() + bborder.x()) - (min_fluid_reco.x() - bborder.x())) /
                          cell_width +
                      1;
            uint ny = ((max_fluid_reco.y() + bborder.y()) - (min_fluid_reco.y() - bborder.y())) /
                          cell_width +
                      1;
            uint nz = ((max_fluid_reco.z() + bborder.z()) - (min_fluid_reco.z() - bborder.z())) /
                          cell_width +
                      1;

            learnSPH::densities::compute_fluid_density_surface_reco(
                fluid_densities_for_surface_reco, particles_positions, point_set_id_fluid, ps_fluid,
                cubic_kernel);

            learnSPH::theta_functions::FluidThetaFunction fluidSDF(cubic_kernel, c, cell_width,
                                                                   beta, nx + 1, ny + 1, nz + 1);
            learnSPH::surface::MarchingCubes mcubes(cell_width, nx, ny, nz,
                                                    min_fluid_reco - bborder, epsilon);

            if (surface_reco_method == 0) {
                std::vector<double> level_set((nx + 1) * (ny + 1) * (nz + 1), -c);
                fluidSDF.computeLevelSet(level_set, particles_positions,
                                         fluid_densities_for_surface_reco,
                                         min_fluid_reco - bborder);
                mcubes.get_Isosurface(level_set);
            } else {
                std::unordered_map<uint64_t, double> level_map;
                fluidSDF.computeLevelMap(level_map, particles_positions,
                                         fluid_densities_for_surface_reco,
                                         min_fluid_reco - bborder);
                mcubes.get_Isosurface_sparse(level_map);
            }

            mcubes.compute_normals();
            write_tri_mesh_to_vtk(filename, mcubes.intersections, mcubes.triangles,
                                  mcubes.intersectionNormals);

            t_next_frame += sim_setup.t_between_frames;

            auto progress_msg = utils::updateProgressBar(stepCounter, maxSteps, 75);
            std::cout << progress_msg.str() << "\n";
            utils::logMessage(progress_msg.str(), log_file);
        }
    }
    return 0;
}