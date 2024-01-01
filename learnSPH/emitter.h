#ifndef EMITTER
#define EMITTER
#include "../extern/CompactNSearch/include/CompactNSearch/CompactNSearch"
#include <Eigen/Dense>
#include <array>
#include <vector>

namespace learnSPH
{
namespace emitter
{
class Emitter
{
  public:
    // emitting direction
    Eigen::Vector3d dir, origin;
    Eigen::Matrix3d rotation_matrix;
    double particle_radius, r, emit_velocity, last_emit, particle_diameter;
    bool alt_pattern, continous;
    int batch_size;
    int emit_counter;
    double r_squared;
    std::vector<std::array<int, 3>> &emit_mark;
    std::vector<Eigen::Vector3d> &particles_positions, &particles_accelerations,
        &particles_velocities;
    std::vector<double> &particles_densities, &particles_pressure;
    std::vector<bool> &deleteFlag;
    CompactNSearch::NeighborhoodSearch &nsearch;
    unsigned int &point_set_id_fluid;

    Emitter(const Eigen::Vector3d dir, const Eigen::Vector3d origin, const double r,
            const double particles_radius, const double emit_velocity, const int emit_counter,
            std::vector<std::array<int, 3>> &emit_mark,
            std::vector<Eigen::Vector3d> &particles_positions,
            std::vector<Eigen::Vector3d> &particles_accelerations,
            std::vector<Eigen::Vector3d> &particles_velocities,
            std::vector<double> &particles_densities, std::vector<double> &particles_pressure,
            std::vector<bool> &deleteFlag, unsigned int &point_set_id_fluid,
            CompactNSearch::NeighborhoodSearch &nsearch);

    void emit_particles(double t_sim, int idx);
    void emit_particles_alternating(double t_sim, int idx);
};
} // namespace emitter
} // namespace learnSPH
#endif