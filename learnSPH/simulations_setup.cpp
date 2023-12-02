#include "simulations_setup.h"
#include <iostream>

learnSPH::simulations_setup::Simulations::Simulations(){}


void learnSPH::simulations_setup::Simulations::simple_cube(){
    this->particle_radius = 0.03;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(1);
    this->fluid_end.resize(1);
    this->fluid_velocities.resize(1);

    this->fluid_begin[0] = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->fluid_end[0] = Eigen::Vector3d(1.0, 1.0, 1.0);
    this->fluid_velocities[0] = Eigen::Vector3d(0.0, 0.0, 0.0);

    this->dt_default = 0.001;
    this->t_between_frames = 0.005;
    this->B = 1000 * 1.02;
    this->v_f = 0;
    this->v_b = 0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->assignment = "assignment2/simple_cube";
}

void learnSPH::simulations_setup::Simulations::simple_cube_with_fluid_viscosity(){
    this->particle_radius = 0.03;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(1);
    this->fluid_end.resize(1);
    this->fluid_velocities.resize(1);

    this->fluid_begin[0] = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->fluid_end[0] = Eigen::Vector3d(1.0, 1.0, 1.0);
    this->fluid_velocities[0] = Eigen::Vector3d(0.0, 0.0, 0.0);

    this->dt_default = 0.001;
    this->t_between_frames = 0.005;
    this->B = 1000 * 1.02;
    this->v_f = 0.1;
    this->v_b = 0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->assignment = "assignment2/simple_cube_with_fluid_viscosity";
}

void learnSPH::simulations_setup::Simulations::cubes_colision(){
    this->particle_radius = 0.03;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(2);
    this->fluid_end.resize(2);
    this->fluid_velocities.resize(2);

    this->fluid_begin[0] = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->fluid_end[0] = Eigen::Vector3d(1.0, 1.0, 1.0);
    this->fluid_velocities[0] = Eigen::Vector3d(2.0, 0.0, 0.0);

    this->fluid_begin[1] = Eigen::Vector3d(1.1, 0.0, 0.0);
    this->fluid_end[1] = Eigen::Vector3d(2.1, 1.0, 1.0);
    this->fluid_velocities[1] = Eigen::Vector3d(-2.0, 0.0, 0.0);

    this->dt_default = 0.001;
    this->t_between_frames = 0.005;
    this->B = 1000 * 1.02;
    this->v_f = 0.1;
    this->v_b = 0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->assignment = "assignment2/cubes_colision";

}

void learnSPH::simulations_setup::Simulations::just_gravity(){
    this->particle_radius = 0.03;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(1);
    this->fluid_end.resize(1);
    this->fluid_velocities.resize(1);

    this->fluid_begin[0] = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->fluid_end[0] = Eigen::Vector3d(1.0, 1.0, 1.0);
    this->fluid_velocities[0] = Eigen::Vector3d(0.0, 0.0, 0.0);

    this->dt_default = 0.001;
    this->t_between_frames = 0.005;
    this->B = 1000 * 1.02;
    this->v_f = 0.1;
    this->v_b = 0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.8);
    this->assignment = "assignment2/just_gravity";
}
void learnSPH::simulations_setup::Simulations::gravity_with_floor(){
    this->particle_radius = 0.03;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(1);
    this->fluid_end.resize(1);
    this->fluid_velocities.resize(1);

    this->fluid_begin[0] = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->fluid_end[0] = Eigen::Vector3d(1.0, 1.0, 1.0);
    this->fluid_velocities[0] = Eigen::Vector3d(0.0, 0.0, 0.0);

    this->boundaries.resize(1);
    this->boundaries[0].filename = "./res/floor.obj";
    this->boundaries[0].min = Eigen::Vector3d(-2, -2, -0.05);
    this->boundaries[0].max = Eigen::Vector3d(2, 2, 0.05);
    this->boundaries[0].inner = false;

    this->dt_default = 0.001;
    this->t_between_frames = 0.005;
    this->B = 1000 * 1.02;
    this->v_f = 0.1;
    this->v_b = 0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.8);
    this->assignment = "assignment2/gravity_with_floor";
}

void learnSPH::simulations_setup::Simulations::gravity_with_floor_boundary_viscosity(){
    this->particle_radius = 0.03;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(1);
    this->fluid_end.resize(1);
    this->fluid_velocities.resize(1);

    this->fluid_begin[0] = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->fluid_end[0] = Eigen::Vector3d(1.0, 1.0, 1.0);
    this->fluid_velocities[0] = Eigen::Vector3d(0.0, 0.0, 0.0);

    this->boundaries.resize(1);
    this->boundaries[0].filename = "./res/floor.obj";
    this->boundaries[0].min = Eigen::Vector3d(-2, -2, -0.05);
    this->boundaries[0].max = Eigen::Vector3d(2, 2, 0.05);
    this->boundaries[0].inner = false;

    this->dt_default = 0.001;
    this->t_between_frames = 0.005;
    this->B = 1000 * 1.02;
    this->v_f = 0.1;
    this->v_b = 1;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.8);
    this->assignment = "assignment2/gravity_with_floor_boundary_viscosity";
}

void learnSPH::simulations_setup::Simulations::dam_break()
{
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;


    this->fluid_begin.resize(1);
    this->fluid_end.resize(1);
    this->fluid_velocities.resize(1);

    this->fluid_begin[0] = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->fluid_end[0] = Eigen::Vector3d(0.15, 0.25, 0.5);
    this->fluid_velocities[0] = Eigen::Vector3d(0.0, 0.0, 0.0);

    this->boundaries.resize(1);
    this->boundaries[0].filename = "./res/boundary.obj";
    this->boundaries[0].min = Eigen::Vector3d(-0.02, -0.02, -0.02);
    this->boundaries[0].max = Eigen::Vector3d(0.17, 0.8, 1.0);
    this->boundaries[0].inner = false;


    this->dt_default = 0.00025;
    this->t_between_frames =  0.006;
    this->B = 1000 * 1.02;
    this->v_f = 0.0025;
    this->v_b = 0.0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
    this->assignment = "assignment2/dam_break";
}

void learnSPH::simulations_setup::Simulations::our_simulation_scene(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(1);
    this->fluid_end.resize(1);
    this->fluid_velocities.resize(1);
    this->fluid_begin[0] = Eigen::Vector3d(0.0, 0.2, 0.4);
    this->fluid_end[0] = Eigen::Vector3d(0.15, 0.6, 0.8);
    this->fluid_velocities[0] = Eigen::Vector3d(0.0, 0.0, 0.0);

    this->boundaries.resize(2);
    this->boundaries[0].filename = "./res/boundary.obj";
    this->boundaries[0].min = Eigen::Vector3d(-0.02, -0.02, -0.02);
    this->boundaries[0].max = Eigen::Vector3d(0.17, 0.8, 1.0);
    this->boundaries[0].inner = false;

    this->boundaries[1].filename = "./res/inner_box.obj";
    this->boundaries[1].min = Eigen::Vector3d(0.05, 0.3, 0.2);
    this->boundaries[1].max = Eigen::Vector3d(0.1, 0.5, 0.3);
    this->boundaries[1].inner = true;

    this->dt_default = 0.00025;
    this->t_between_frames = 0.006;
    this->B = 1000 * 1.02;
    this->v_f = 0.0025;
    this->v_b = 0.0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
    this->assignment = "assignment3/our_simulation";
}