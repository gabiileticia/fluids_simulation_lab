#include "acceleration.h"
#include <stdlib.h> // rand
#include <omp.h>
#include "kernel.h"
#include <algorithm>
#include <vector>

learnSPH::acceleration::Acceleration::Acceleration(
        double B
        , double v_f
        , double v_b
        , double h
        , double roh_0
        , Eigen::Vector3d gravity
        , learnSPH::kernel::CubicSplineKernel &kernel) : B(B)
                                                      , v_f(v_f)
                                                      , v_b(v_b)
                                                      , h(h)
                                                      , kernel(kernel)
                                                      , roh_0(roh_0)
                                                      , gravity(gravity)
{
    this->h_square = h * h;
}

void learnSPH::acceleration::Acceleration::pressure(
        std::vector<double> &particles_pressure
        , std::vector<double> &particles_density
        , double rest_density) 
{
    double delta_density;

    for (int i = 0; i < particles_density.size(); i++) {
        delta_density = this->B * (particles_density[i] - rest_density);
        particles_pressure[i] = std::max(0.0, delta_density);
    }
}

void learnSPH::acceleration::Acceleration::accelerations(
        std::vector<Eigen::Vector3d> &accelerations,
        std::vector<double> &roh, // density
        std::vector<double> &p,     // pressure
        unsigned int ps_id_fluid,
        unsigned int ps_id_boundary,
        CompactNSearch::PointSet const& fluid_neighbors,
        CompactNSearch::PointSet const& boundary_neighbors,
        std::vector<Eigen::Vector3d> &fluid_particles,
        std::vector<Eigen::Vector3d> &boundary_particles,
        std::vector<Eigen::Vector3d> &velocity,
        double boundary_volume,
        const double fluid_mass)
{

    Eigen::Vector3d ap, av, ae; // pressure, velocity and mass forces
    Eigen::Vector3d ff_inter_pressure
                , fs_inter_pressure
                , ff_inter_velocity
                , fs_inter_velocity; // fluid fluid and fluid solid interaction
    Eigen::Vector3d dx;
    Eigen::Vector3d kernel_grad;
    double dx_norm;
    double h_square_cent = 0.01 * this->h_square;
    double roh_square_i_inverse;

    for (int i = 0; i < fluid_neighbors.n_points(); ++i) {

        ff_inter_pressure = {0, 0, 0};
        fs_inter_pressure = {0, 0, 0};
        ff_inter_velocity = {0, 0, 0};
        fs_inter_velocity = {0, 0, 0};
        roh_square_i_inverse = 1 / (roh[i] * roh[i]);
        
        for (size_t j = 0; j < fluid_neighbors.n_neighbors(ps_id_fluid, i); j++) {
  
            unsigned int pid = fluid_neighbors.neighbor(ps_id_fluid, i, j);
            dx = fluid_particles[i] - fluid_particles[pid];
            kernel_grad = kernel.kernel_gradient(dx);
            dx_norm = dx.norm();

            ff_inter_pressure += fluid_mass *
                    (p[i] * roh_square_i_inverse + p[pid] / (roh[pid] * roh[pid])) * kernel_grad;            

            ff_inter_velocity += fluid_mass * 
                    (velocity[i] - velocity[pid]) * (dx).transpose() * kernel_grad / (roh[pid] * (dx_norm * dx_norm + h_square_cent));
        }

        
        for (size_t j = 0; j < fluid_neighbors.n_neighbors(ps_id_boundary, i); ++j) 
        {
            unsigned int pid = fluid_neighbors.neighbor(ps_id_boundary, i, j);
            dx = fluid_particles[i] - boundary_particles[pid];
            kernel_grad = kernel.kernel_gradient(dx);
            dx_norm = dx.norm();

            fs_inter_pressure += this-> roh_0 * boundary_volume * p[i] * roh_square_i_inverse * kernel_grad;

            fs_inter_velocity += boundary_volume * velocity[i] * dx.transpose() * kernel_grad / (dx_norm * dx_norm + h_square_cent);
        }

        ap = -ff_inter_pressure - fs_inter_pressure;
        av = 2 * this->v_f * ff_inter_velocity + 2 * this->v_b * fs_inter_velocity;
        ae = this->gravity;

        accelerations[i] = ap + av + ae;
    }
}