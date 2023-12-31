#include "emitter.h"
#include "utils.h"
#include <array>
#include <vector>

learnSPH::emitter::Emitter::Emitter(
    const Eigen::Vector3d dir, const Eigen::Vector3d origin, const double r,
    const double particle_radius, const double emit_velocity,
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
    // Calculate rotation matrix to align direction with positive z-axis
    Eigen::Vector3d unit_z(0, 0, 1); // Positive z-axis
    Eigen::Vector3d axis = this->dir.cross(unit_z);
    float angle          = std::acos(this->dir.dot(unit_z) / (this->dir.norm() * unit_z.norm()));
    Eigen::AngleAxisd rotation(angle, axis.normalized());
    Eigen::Matrix3d rotation_matrix = rotation.toRotationMatrix();

    // calculation translation to new origin
    Eigen::Matrix4d translation_matrix   = Eigen::Matrix4d::Identity();
    translation_matrix.block<3, 1>(0, 3) = this->origin;

    // calculate final transformation matrix
    transformation_matrix                   = Eigen::Matrix4d::Identity();
    transformation_matrix.block<3, 3>(0, 0) = rotation_matrix;
    transformation_matrix *= translation_matrix;

    this->particle_diameter = this->particle_radius * 2;
    this->last_emit         = 0;
}

void learnSPH::emitter::Emitter::emit_particles(double t_sim, int idx)
{
    double corner;
    double x, y;
    Eigen::Vector4d new_particle;
    std::array<int, 3> new_batch;

    this->last_emit = t_sim;

    int new_part_counter = 0;
    int l                = (int)((2 * this->r) / this->particle_diameter);
    std::vector<Eigen::Vector3d> new_particles;
    new_particles.resize(l * l);

    x      = 0;
    y      = 0;
    corner = -l * this->particle_diameter;

    for (int i = 0; i < l; i++) {
        x = corner + this->particle_diameter * i;
        for (int j = 0; j < l; j++) {
            y = corner + this->particle_diameter * j;
            if (((x * x + y * y) - (this->r)) < 0) {
                new_particle                    = {x, y, 0, 1};
                new_particle                    = transformation_matrix * new_particle;
                new_particles[new_part_counter] = new_particle.head<3>();
                new_part_counter++;
            }
        }
        y = 0;
    }
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
        this->particles_velocities[new_batch[0] + i] = this->dir * this->emit_velocity;
    }
}

void learnSPH::emitter::Emitter::emit_particles_alternating(double t_sim, int idx){};