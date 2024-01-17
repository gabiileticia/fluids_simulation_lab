#include "emitter.h"
#include "Eigen/src/Core/Matrix.h"
#include "Eigen/src/Geometry/AngleAxis.h"
#include "utils.h"
#include <array>
#include <cmath>
#include <vector>

learnSPH::emitter::Emitter::Emitter(
    const Eigen::Vector3d dir, const Eigen::Vector3d origin, const double r,
    const double particle_radius, const double emit_velocity, const int emit_counter,
    std::vector<std::array<int, 3>> &emit_mark, std::vector<Eigen::Vector3d> &particles_positions,
    std::vector<Eigen::Vector3d> &particles_accelerations,
    std::vector<Eigen::Vector3d> &particles_velocities, std::vector<double> &particles_densities,
    std::vector<double> &particles_pressure, std::vector<bool> &deleteFlag,
    unsigned int &point_set_id_fluid, CompactNSearch::NeighborhoodSearch &nsearch)
    : dir(dir), origin(origin), r(r), particle_radius(particle_radius),
      emit_velocity(emit_velocity), emit_mark(emit_mark), particles_positions(particles_positions),
      particles_accelerations(particles_accelerations), particles_velocities(particles_velocities),
      particles_densities(particles_densities), particles_pressure(particles_pressure),
      deleteFlag(deleteFlag), point_set_id_fluid(point_set_id_fluid), nsearch(nsearch)
{
    this->dir.normalize();
    Eigen::Vector3d unit_z = {0, 0, 1};
    unit_z.normalize();

    Eigen::Vector3d axis = unit_z.cross(this->dir);
    double angle         = std::acos(unit_z.dot(this->dir));
    axis.normalize();

    Eigen::AngleAxisd rotation(angle, axis);
    rotation_matrix = rotation.toRotationMatrix();
    this->alt_pattern = false;

    this->particle_diameter = this->particle_radius * 2;
    this->last_emit         = 0;
    this->r_squared         = this->r * this->r;
    if (emit_counter == 0) {
        this->continous    = true;
        this->emit_counter = 1;
    } else {
        this->continous    = false;
        this->emit_counter = emit_counter;
    }
}

void learnSPH::emitter::Emitter::emit_particles(double t_sim, int idx)
{
    double corner;
    double x, y;
    Eigen::Vector3d new_particle;
    std::array<int, 3> new_batch;

    this->last_emit = t_sim;

    int new_part_counter = 0;
    int l                = (int)(this->r/ this->particle_diameter);
    std::vector<Eigen::Vector3d> new_particles;
    double new_r_squared = l * this->particle_radius * l * this->particle_radius;
    new_particles.resize(l * l * 4);

    double sqrt3 = std::sqrt(3);
    double onethird = 1. / 3.;

    x      = 0;
    y      = 0;
    corner = -l * this->particle_diameter;

    for (int i = 0; i < 2*l; i++) {
        x = corner + this->particle_diameter * i;
        for (int j = 0; j < 2*l; j++) {
            y = corner + this->particle_diameter * j;
            if (((x * x + y * y) - new_r_squared )< 0) {
                new_particle = this->rotation_matrix * Eigen::Vector3d({x, y, 0}) + this->origin;
                new_particles[new_part_counter] = new_particle;
                new_part_counter++;
            }
        }
        y = 0;
    }
    new_particles.resize(new_part_counter);
    // std::cout << "pos before " << this->particles_positions.size() << "\n";
    new_batch[0] = this->particles_positions.size();
    new_batch[1] = this->particles_positions.size() + new_part_counter;
    new_batch[2] = idx;
    this->emit_mark.push_back(new_batch);

    this->particles_positions.resize(this->particles_positions.size() + new_particles.size());
    // std::cout << "pos after " << this->particles_positions.size() << "\n";
    this->particles_accelerations.resize(this->particles_accelerations.size() +
                                         new_particles.size());
    this->particles_velocities.resize(this->particles_velocities.size() + new_particles.size());
    this->particles_densities.resize(this->particles_densities.size() + new_particles.size());
    this->particles_pressure.resize(this->particles_pressure.size() + new_particles.size());
    this->deleteFlag.resize(this->deleteFlag.size() + new_particles.size());

    this->nsearch.resize_point_set(this->point_set_id_fluid,
                                   this->particles_positions.front().data(),
                                   this->particles_positions.size());

    for (int i = 0; i < new_part_counter; i++) {
        this->particles_positions[new_batch[0] + i]  = new_particles[i];
        this->particles_velocities[new_batch[0] + i] = this->dir.normalized() * this->emit_velocity;
        //this->particles_accelerations[new_batch[0] + i] = this->dir.normalized() * this->emit_velocity / 2;
    }
    if (!(this->continous))
        this->emit_counter--;
}

void learnSPH::emitter::Emitter::emit_particles_alternating(double t_sim, int idx){
    double corner;
    double x, y;
    Eigen::Vector3d new_particle;
    std::array<int, 3> new_batch;

    this->last_emit = t_sim;

    int new_part_counter = 0;
    int l                = (int)(this->r/ this->particle_diameter);
    std::vector<Eigen::Vector3d> new_particles;
    double new_r_squared = l * this->particle_radius * l * this->particle_radius;
    new_particles.resize(l * l * 4);

    double sqrt3 = std::sqrt(3);
    double onethird = 1. / 3.;

    x      = 0;
    y      = 0;
    corner = -l * this->particle_diameter;

    for (int i = 0; i < 2*l; i++) {
        for (int j = 0; j < 2*l; j++) {
            x = corner + this->particle_radius * (2*i + ((j + this->alt_pattern)%2));
            y = corner + this->particle_radius * (sqrt3 * (j + onethird * this->alt_pattern));
            if (((x * x + y * y) - new_r_squared )< 0) {
                new_particle = this->rotation_matrix * Eigen::Vector3d({x, y, 0}) + this->origin;
                new_particles[new_part_counter] = new_particle;
                new_part_counter++;
            }
        }
        y = 0;
    }
    // for AB AB AB AB pattern in HCP generation
    this->alt_pattern ^= true;

    new_particles.resize(new_part_counter);
    // std::cout << "pos before " << this->particles_positions.size() << "\n";
    new_batch[0] = this->particles_positions.size();
    new_batch[1] = this->particles_positions.size() + new_part_counter;
    new_batch[2] = idx;
    this->emit_mark.push_back(new_batch);

    this->particles_positions.resize(this->particles_positions.size() + new_particles.size());
    // std::cout << "pos after " << this->particles_positions.size() << "\n";
    this->particles_accelerations.resize(this->particles_accelerations.size() +
                                         new_particles.size());
    this->particles_velocities.resize(this->particles_velocities.size() + new_particles.size());
    this->particles_densities.resize(this->particles_densities.size() + new_particles.size());
    this->particles_pressure.resize(this->particles_pressure.size() + new_particles.size());
    this->deleteFlag.resize(this->deleteFlag.size() + new_particles.size());

    this->nsearch.resize_point_set(this->point_set_id_fluid,
                                   this->particles_positions.front().data(),
                                   this->particles_positions.size());

    for (int i = 0; i < new_part_counter; i++) {
        this->particles_positions[new_batch[0] + i]  = new_particles[i];
        this->particles_velocities[new_batch[0] + i] = this->dir.normalized() * this->emit_velocity;
        //this->particles_accelerations[new_batch[0] + i] = this->dir.normalized() * this->emit_velocity / 2;
    }
    if (!(this->continous))
        this->emit_counter--;
};