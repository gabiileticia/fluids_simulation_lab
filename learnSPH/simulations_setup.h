#include <Eigen/Dense>
#include <vector>


#ifndef SIMULATION
#define SIMULATION

namespace learnSPH {
    namespace simulations_setup {
        class Simulations {
            public:
            double particle_radius;

            double fluid_rest_density;

            std::vector<Eigen::Vector3d> fluid_begin;
            std::vector<Eigen::Vector3d> fluid_end;
            std::vector<Eigen::Vector3d> fluid_velocities;

            std::vector<std::string> boundary_file;
            Eigen::Vector3d boundary_begin;
            Eigen::Vector3d boundary_end;

            double dt_default;
            double t_between_frames;
            double B;
            double v_f;
            double v_b;
            Eigen::Vector3d gravity;
            std::string assignment;

            Simulations(std::string teste);

            void simple_cube();
            void simple_cube_with_fluid_viscosity();
            void cubes_colision();
            void just_gravity();
            void gravity_with_floor();
            void gravity_with_floor_boundary_viscosity();
            void dam_break();
            void our_simulation_scene();

        };
    }
}

#endif





