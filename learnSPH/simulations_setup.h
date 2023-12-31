// #include <Eigen/Dense>
#include <vector>

#include <Eigen/Dense>
#include "types.h"
#include "../learnSPH/emitter.h"


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
            std::vector<types::object> objects;

            std::vector<learnSPH::types::emitter_data> emitters;

            double dt_default;
            double t_between_frames;
            double B;
            double v_f;
            double v_b;
            Eigen::Vector3d gravity;
            std::string assignment;
            bool simbound_active;
            Eigen::Vector3d sim_boundary_min;
            Eigen::Vector3d sim_boundary_max;

            bool no_fluid = false;

            Simulations();

            void simple_cube();
            void simple_cube_with_fluid_viscosity();
            void cubes_colision();
            void just_gravity();
            void gravity_with_floor();
            void gravity_with_floor_boundary_viscosity();
            void dam_break();
            void our_simulation_scene();
            void mcubes_stress_test_scene();
            void slope_ramp_wall_vessel();
            void empty_scene_test();
            void simple_emitter_test();

        };
    }
}

#endif






