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
#include "../learnSPH/utils.h"

#include <chrono>

int main() {
  std::cout << "Welcome to the learnSPH framework!!" << std::endl;
  std::cout << "Yor current setup is:" << std::endl;

  // Initializing
  double particle_radius = 0.005;
  std::cout << "particle_radius: " << particle_radius << std::endl;
  const double particle_diameter = 2.0 * particle_radius;
  std::cout << "particle_diameter: " << particle_diameter << std::endl;
  const double fluid_sampling_distance = particle_diameter;
  std::cout << "fluid_sampling_distance: " << fluid_sampling_distance
            << std::endl;
  const double boundary_sampling_distance = 0.8 * particle_diameter;
  std::cout << "boundary_sampling_distance: " << boundary_sampling_distance
            << std::endl;

  // smoothing lenght
  const double h = 1.2 * particle_diameter;
  std::cout << "h: " << h << std::endl;
  // compact support
  const double beta = 2.0 * h;
  std::cout << "beta: " << beta << std::endl;

  const double fluid_rest_density = 1000.0;
  std::cout << "fluid_density: " << fluid_rest_density << std::endl;
  Eigen::Vector3d fluid_begin = Eigen::Vector3d(0.01, 0.01, 0.01);
  std::cout << "fluid_begin: " << fluid_begin.transpose() << std::endl;
  Eigen::Vector3d fluid_end = Eigen::Vector3d(0.125, 0.25, 0.5);
  std::cout << "fluid_end: " << fluid_end.transpose() << std::endl;
  const double fluid_volume = fluid_end.x() * fluid_end.y() * fluid_end.z();
  std::cout << "fluid_volume: " << fluid_volume << std::endl;
  const double fluid_mass = fluid_volume * fluid_rest_density;
  std::cout << "fluid_mass: " << fluid_mass << std::endl;

  const std::string boundary_file = "./res/boundary.obj";
  std::cout << "boundary_file: " << boundary_file << std::endl;
  Eigen::Vector3d boundary_begin = Eigen::Vector3d(0.0, 0.0, 0.0);
  std::cout << "boundary_begin: " << boundary_begin.transpose() << std::endl;
  Eigen::Vector3d boundary_end = Eigen::Vector3d(0.15, 0.8, 1.0);
  std::cout << "boundary_end: " << boundary_end.transpose() << std::endl;
  double a, b, c;
  a = 0.2;
  b = 1.0;
  c = 1.0;
  // const double particle_volume = 4.0/3.0 * learnSPH::kernel::PI *
  // particle_radius*particle_radius*particle_radius; const double
  // boundary_volume = 2 * particle_volume * (a*b + a*c + a*c);
  const double boundary_volume =
      0.001 * boundary_end.x() * boundary_end.y() * boundary_end.z();
  // const double boundary_volume = particle_diameter * (2 * boundary_end.x() *
  // boundary_end.y() + 													2 * boundary_end.x() * boundary_end.z() + 													2 *
  // boundary_end.y() * boundary_end.z());
  std::cout << "boundary_volume: " << boundary_volume << std::endl;
  const double boundary_mass = boundary_volume * fluid_rest_density;
  std::cout << "boundary_mass: " << boundary_mass << std::endl;

  double dt_cfl, dt;
  double dt_default = 0.00025;
  std::cout << "dt_default: " << dt_default << std::endl;
  double t_next_frame = 0;
  double t_between_frames = 0.0005;
  double t_simulation = 0;
  double B = 1000 * 1.02;
  std::cout << "B: " << B << std::endl;
  double v_f = 0.0025;
  std::cout << "v_f: " << v_f << std::endl;
  double v_b = 0.0;
  std::cout << "v_b: " << v_b << std::endl;
  Eigen::Vector3d gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
  std::cout << "gravity: " << gravity.transpose() << std::endl;

  // instantiating some classes
  learnSPH::kernel::CubicSplineKernel cubic_kernel(h);
  learnSPH::acceleration::Acceleration acceleration(
      B, v_f, v_b, h, fluid_rest_density, gravity, cubic_kernel);
  learnSPH::timeIntegration::semiImplicitEuler semImpEuler(
      particle_radius, boundary_end, boundary_begin);

  // Load simulation geometry
  std::vector<Eigen::Vector3d> boundary_particles_positions;
  learnSPH::geometry::load_n_sample_boundary(
      boundary_particles_positions, boundary_file, boundary_sampling_distance);

  std::cout << "Number of boundary particles" << std::endl;
  std::cout << boundary_particles_positions.size() << std::endl;

  // Sample the fluid with particles
  std::vector<Eigen::Vector3d> particles_positions;
  learnSPH::sampling::fluid_box(particles_positions, fluid_begin, fluid_end,
                                fluid_sampling_distance);

  // int end_of_box1 = particles_positions.size() - 1;

  // learnSPH::sampling::fluid_box(particles_positions, {1.01, 0.0, 0.0},
  //                             {2.01, 1.0, 1.0}, fluid_sampling_distance);
  // int end_of_box2 = particles_positions.size() - 1;

  std::cout << "Number of fluid particles" << std::endl;
  std::cout << particles_positions.size() << std::endl;

  double fluid_particle_mass = fluid_mass / particles_positions.size();
  std::cout << "fluid_particle_mass: " << fluid_particle_mass << std::endl;

  std::vector<double> particles_densities(particles_positions.size());
  std::vector<Eigen::Vector3d> particles_accelerations(
      particles_positions.size());
  std::vector<double> particles_pressure(particles_positions.size());

  std::vector<Eigen::Vector3d> particles_velocities(
      particles_positions.size(), Eigen::Vector3d(0.0, 0.0, 0.0));

  // for (int i = 0; i <= end_of_box1; ++i) {
  //     particles_velocities[i] = {2, 0, 0};
  // }
  // for (int i = end_of_box1 + 1; i <= end_of_box2; ++i) {
  //     particles_velocities[i] = {-2, 0, 0};
  // }

  // Setting up neighborhood search
  CompactNSearch::NeighborhoodSearch nsearch(beta);
  unsigned int point_set_id_boundary =
      nsearch.add_point_set(boundary_particles_positions.front().data(),
                            boundary_particles_positions.size());

  unsigned int point_set_id_fluid = nsearch.add_point_set(
      particles_positions.front().data(), particles_positions.size());

  nsearch.set_active(point_set_id_boundary, point_set_id_fluid, false);
  nsearch.set_active(point_set_id_boundary, point_set_id_boundary, true);
  nsearch.set_active(point_set_id_fluid, point_set_id_fluid, true);
  nsearch.set_active(point_set_id_fluid, point_set_id_boundary, true);
  nsearch.find_neighbors();

  CompactNSearch::PointSet const &ps_boundary =
      nsearch.point_set(point_set_id_boundary);
  CompactNSearch::PointSet const &ps_fluid =
      nsearch.point_set(point_set_id_fluid);

  // Compute boundary masses
  const double boundary_particles_mass =
      boundary_mass / boundary_particles_positions.size();
  const double boundary_particle_density =
      fluid_rest_density / particles_positions.size();
  std::cout << "boundary_particles_mass: " << boundary_particles_mass
            << std::endl;

  std::vector<double> boundary_particles_densities(
      boundary_particles_positions.size());
  learnSPH::densities::compute_boundary_densities(
      boundary_particles_densities, boundary_particles_positions,
      point_set_id_boundary, ps_boundary, boundary_particles_mass,
      cubic_kernel);

  std::vector<double> boundary_particles_masses(
      boundary_particles_positions.size());
  learnSPH::densities::compute_boundary_masses(
      boundary_particles_masses, boundary_particles_positions,
      point_set_id_boundary, ps_boundary, fluid_rest_density, cubic_kernel);

  int total_vector =
      particles_positions.size() + boundary_particles_positions.size();
  std::vector<Eigen::Vector3d> output_positions(total_vector);
  std::vector<double> output_densities(total_vector);
  std::vector<double> b_densities(boundary_particles_positions.size(),
                                  fluid_rest_density);
  std::vector<Eigen::Vector3d> output_velocities(total_vector);
  std::vector<Eigen::Vector3d> b_velocities(boundary_particles_positions.size(),
                                            Eigen::Vector3d(0.0, 0.0, 0.0));

  // keeping track of number of elements which will be deleted
  int count_del = 0;
  std::vector<bool> deleteFlag(particles_positions.size());

  // Simulation loop
  while (t_simulation < 5) {

    std::cout << "t_simulation " << t_simulation << std::endl;
    std::cout << "v_max " << semImpEuler.v_max << std::endl;
    // Compute dt
    dt_cfl = 0.5 * particle_radius * (1 / std::min(50.0, semImpEuler.v_max));
    dt = std::min(dt_cfl, dt_default);
    // std::cout << "dt " << dt << std::endl;

    // Find neighbors
    nsearch.find_neighbors();

    // Compute fluid particles densities
    learnSPH::densities::compute_fluid_density(
        particles_densities, particles_positions, boundary_particles_positions,
        boundary_particles_masses, point_set_id_fluid, ps_fluid,
        point_set_id_boundary, ps_boundary, fluid_particle_mass, cubic_kernel);

    // Compute acceleration
    acceleration.pressure(particles_pressure, particles_densities,
                          fluid_rest_density);

    acceleration.accelerations(
        particles_accelerations, particles_densities, particles_pressure,
        point_set_id_fluid, point_set_id_boundary, ps_fluid, ps_boundary,
        particles_positions, boundary_particles_positions, particles_velocities,
        boundary_particles_masses, boundary_particle_density,
        fluid_particle_mass);

    // Integrate
    semImpEuler.integrationStep(particles_positions, particles_velocities,
                                particles_accelerations, deleteFlag, dt,
                                count_del);
    learnSPH::utils::deleteOutOfBounds(
        particles_positions, particles_velocities, particles_accelerations,
        particles_densities, particles_pressure, deleteFlag, count_del);
    nsearch.resize_point_set(point_set_id_fluid, ps_fluid, particles_positions.size());

    auto stop_integrate = std::chrono::high_resolution_clock::now();

    // Increment t
    t_simulation += dt;

    // Save output
    if (t_simulation >= t_next_frame) {
      const std::string filename =
          "./res/wcsph_dam_break/wcsph_" +
          std::to_string((int)(t_simulation * 1000000)) + ".vtk";

      learnSPH::write_particles_to_vtk(filename, particles_positions,
                                       particles_pressure,
                                       particles_velocities);
      t_next_frame += t_between_frames;
    }
  }
  return 0;
}