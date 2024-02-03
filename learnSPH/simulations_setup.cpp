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

    this->objects.resize(1);
    this->objects[0].filename = "./res/floor.obj";
    this->objects[0].min = Eigen::Vector3d(-2, -2, -0.05);
    this->objects[0].max = Eigen::Vector3d(2, 2, 0.05);
    this->objects[0].noCheck = true;

    this->sim_boundary_min = Eigen::Vector3d(-2, -2, -0.05);
    this->sim_boundary_max = Eigen::Vector3d(2, 2, 1.05);
    this->simbound_active = true;

    this->dt_default = 0.001;
    this->t_between_frames = 0.005;
    this->B = 1000 * 1.02;
    this->v_f = 0.1;
    this->v_b = 0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.8);
    this->assignment = "assignment2/gravity_with_floor";

    this->pressure_solver_method = 0;
    this->surface_reco_method = 1;
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

    this->objects.resize(1);
    this->objects[0].filename = "./res/floor.obj";
    this->objects[0].min = Eigen::Vector3d(-2, -2, -0.05);
    this->objects[0].max = Eigen::Vector3d(2, 2, 0.05);
    this->objects[0].noCheck = true;

    this->sim_boundary_min = this->objects[0].min;
    this->sim_boundary_max = this->objects[0].max;
    this->simbound_active = true;

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

    this->objects.resize(1);
    this->objects[0].filename = "./res/boundary.obj";
    this->objects[0].min = Eigen::Vector3d(-0.02, -0.02, -0.02);
    this->objects[0].max = Eigen::Vector3d(0.17, 0.8, 1.0);
    this->objects[0].noCheck = true;

    this->sim_boundary_min = this->objects[0].min;
    this->sim_boundary_max = this->objects[0].max;
    this->simbound_active = true;

    this->surface_reco_method = 1;  // 0: dense; 1: sparse
    this->pressure_solver_method = 1;   // 0: wcsph, 1: pbf
    this->n_iterations_pbf = 5;

    this->dt_default = 0.0002;
    this->t_between_frames =  0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.0025;
    this->v_b = 0.0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
    this->assignment = "assignment4/dam_break_pbf";
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

    this->objects.resize(2);
    this->objects[0].filename = "./res/boundary.obj";
    this->objects[0].min = Eigen::Vector3d(-0.02, -0.02, -0.02);
    this->objects[0].max = Eigen::Vector3d(0.17, 0.8, 1.0);
    this->objects[0].noCheck = true;

    this->sim_boundary_min = this->objects[0].min;
    this->sim_boundary_max = this->objects[0].max;
    this->simbound_active = true;

    this->objects[1].filename = "./res/inner_box.obj";
    this->objects[1].min = Eigen::Vector3d(0.05, 0.3, 0.2);
    this->objects[1].max = Eigen::Vector3d(0.1, 0.5, 0.3);
    this->objects[1].noCheck = false;

    this->dt_default = 0.00025;
    this->t_between_frames = 0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.0025;
    this->v_b = 0.0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
    this->assignment = "our_simulation";
}

void learnSPH::simulations_setup::Simulations::mcubes_stress_test_scene(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(1);
    this->fluid_end.resize(1);
    this->fluid_velocities.resize(1);
    this->fluid_begin[0] = Eigen::Vector3d(0.25,0.25,0.1);
    this->fluid_end[0] = Eigen::Vector3d(0.75,0.75,0.6);
    this->fluid_velocities[0] = Eigen::Vector3d(0,0,0);

    this->objects.resize(1);
    this->objects[0].filename = "./res/boundary_cube.obj";
    this->objects[0].min = Eigen::Vector3d(-0.02,-0.02,-0.02);
    this->objects[0].max = Eigen::Vector3d(1.,1.,1.);
    this->objects[0].noCheck = true;

    this->sim_boundary_min = this->objects[0].min;
    this->sim_boundary_max = this->objects[0].max;
    this->simbound_active = true;

    this->dt_default = 0.00025;
    this->t_between_frames = 0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.0025;
    this->v_b = 0.0;
    this->gravity = Eigen::Vector3d(0.0,0.0,-9.81);
    this->assignment = "mcubes_benchmark";
}

void learnSPH::simulations_setup::Simulations::slope_ramp_wall_vessel(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(1);
    this->fluid_end.resize(1);
    this->fluid_velocities.resize(1);
    this->fluid_begin[0] = Eigen::Vector3d(-0.852252,-0.21327, 1.512171);
    this->fluid_end[0] = Eigen::Vector3d(-0.451804,  0.187177, 1.912618);

    // maybe refactor to give obstacles different name
    this->objects.resize(2);
    this->objects[0].filename = "./res/scene_slope.obj";
    this->objects[0].noCheck = true;
    this->objects[1].filename = "./res/scene_wall_vessel.obj";
    this->objects[1].noCheck = true;

    // simulation domain boundary
    this->sim_boundary_min = Eigen::Vector3d(-1,-1.1,-0.2);
    this->sim_boundary_max = Eigen::Vector3d(2, 1.1, 2.2);
    this->simbound_active = true;

    this->dt_default = 0.001;
    this->t_between_frames = 0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.0025;
    this->v_b = 0.0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
    this->assignment = "assignment4/slope_wall_vessel";
}

void learnSPH::simulations_setup::Simulations::empty_scene_test(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(0);
    this->fluid_end.resize(0);
    this->fluid_velocities.resize(0);

    this->objects.resize(1);
    this->objects[0].filename = "./res/scene_slope.obj";
    this->objects[0].noCheck = true;

    // simulation domain boundary
    this->sim_boundary_min = Eigen::Vector3d(-1,-1.1,-0.2);
    this->sim_boundary_max = Eigen::Vector3d(2, 1.1, 2.2);
    this->simbound_active = true;

    this->dt_default = 0.001;
    this->t_between_frames = 0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.0025;
    this->v_b = 0.0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
    this->assignment = "assignment5/empty_scene_test";
}

void learnSPH::simulations_setup::Simulations::simple_emitter_test(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(0);
    this->fluid_end.resize(0);
    this->fluid_velocities.resize(0);

    this->objects.resize(1);
    this->objects[0].filename = "./res/vessel.obj";
    this->objects[0].noCheck = true;

    // simulation domain boundary
    this->sim_boundary_min = Eigen::Vector3d(0,0,0);
    this->sim_boundary_max = Eigen::Vector3d(1, 1, .5);
    this->simbound_active = true;

    this->surface_reco_method = 1;  // 0: dense; 1: sparse
    this->pressure_solver_method = 1;   // 0: wcsph, 1: pbf
    this->n_iterations_pbf = 10;

    this->dt_default = 0.005;
    this->t_between_frames = 0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.0025;
    this->v_b = 0.0;
    // this->v_f = 0.5;
    // this->v_b = 0.0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
    this->assignment = "assignment5/simple_emitter_test";

    this->emitters.resize(1);
    this->emitters[0].dir = {0,1,0};
    this->emitters[0].origin = {.5,.1,.4};
    this->emitters[0].r = 0.05;
    this->emitters[0].velocity = 1;
    this->emitters[0].alternating = false;
    this->emitters[0].emission_freq = 1;
    this->emitters[0].emit_counter = 1000;
}

void learnSPH::simulations_setup::Simulations::fountain(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(0);
    this->fluid_end.resize(0);
    this->fluid_velocities.resize(0);

    this->objects.resize(1);
    this->objects[0].filename = "./res/vessel.obj";
    this->objects[0].noCheck = true;

    // simulation domain boundary
    this->sim_boundary_min = Eigen::Vector3d(-0.1,-0.1,-0.1);
    this->sim_boundary_max = Eigen::Vector3d(1, 1, 1);
    this->simbound_active = true;

    this->surface_reco_method = 1;  // 0: dense; 1: sparse
    this->pressure_solver_method = 0;   // 0: wcsph, 1: pbf
    this->n_iterations_pbf = 10;

    this->dt_default = 0.00025;
    this->t_between_frames = 0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.0025;
    this->v_b = 0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
    this->assignment = "assignment5/fountain";

    this->emitters.resize(1);
    this->emitters[0].dir = {0,0,1};
    this->emitters[0].origin = {.5,.5,.125};
    this->emitters[0].r = 0.04;
    this->emitters[0].velocity = 2.5;
    this->emitters[0].alternating = true;
    this->emitters[0].emission_freq = 1;
    this->emitters[0].emit_counter = 2000;
}

void learnSPH::simulations_setup::Simulations::water_droplet_no_gravity(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(1);
    this->fluid_end.resize(1);
    this->fluid_velocities.resize(1);

    this->fluid_begin[0] = Eigen::Vector3d(0.5, 0.5, 0.5);
    this->fluid_end[0] = Eigen::Vector3d(0.6, 0.6, 0.6);
    this->fluid_velocities[0] = Eigen::Vector3d(0.0, 0.0, 0.0);

    this->surface_reco_method = 1;  // 0: dense; 1: sparse
    this->pressure_solver_method = 1;   // 0: wcsph, 1: pbf
    this->n_iterations_pbf = 10;

    this->dt_default = 0.0001;
    this->t_between_frames = 0.05;
    this->B = 1000 * 1.02;
    this->v_f = 0.1;
    this->v_b = 0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->assignment = "assignment5/water_droplet";
    this->simTime = 20;

    this->surface_tension = true;
    this->cohesion_coefficient = 0.05;
    this->adhesion_coefficient = 0.01;
}


void learnSPH::simulations_setup::Simulations::boundary_wetting_no_surface_tension(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(0);
    this->fluid_end.resize(0);
    this->fluid_velocities.resize(0);

    this->objects.resize(1);
    this->objects[0].filename = "./res/boundary_cube.obj";
    this->objects[0].min = Eigen::Vector3d(-0.02,-0.02,-0.02);
    this->objects[0].max = Eigen::Vector3d(1.,1.,1.);
    this->objects[0].noCheck = true;

    // this->objects[1].filename = "./res/inner_box.obj";
    // this->objects[1].min = Eigen::Vector3d(0.05, 0.3, 0.2);
    // this->objects[1].max = Eigen::Vector3d(0.1, 0.5, 0.3);
    // this->objects[1].noCheck = false;

    // simulation domain boundary
    this->sim_boundary_min = Eigen::Vector3d(-1,-1,-1);
    this->sim_boundary_max = Eigen::Vector3d(1, 1, 1);
    this->simbound_active = true;

    this->surface_reco_method = 1;  // 0: dense; 1: sparse
    this->pressure_solver_method = 0;   // 0: wcsph, 1: pbf
    this->n_iterations_pbf = 10;

    this->dt_default = 0.0005;
    this->t_between_frames = 0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.1;
    this->v_b = 0.0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
    this->assignment = "assignment5/boundary_wetting";

    this->emitters.resize(1);
    this->emitters[0].dir = {0,1,0};
    this->emitters[0].origin = {.1,0.15,.4};
    this->emitters[0].r = 0.1;
    this->emitters[0].velocity = 0.7;
    this->emitters[0].alternating = false;
    this->emitters[0].emission_freq = 1;
    this->emitters[0].emit_counter = 1000;    


    this->surface_tension = true;
    this->cohesion_coefficient = 0.0; //0.05;
    this->adhesion_coefficient = 0.0; //0.01;
}

void learnSPH::simulations_setup::Simulations::boundary_wetting_only_adhesion(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(0);
    this->fluid_end.resize(0);
    this->fluid_velocities.resize(0);

    this->objects.resize(2);
    this->objects[0].filename = "./res/boundary_cube.obj";
    this->objects[0].min = Eigen::Vector3d(-0.02,-0.02,-0.02);
    this->objects[0].max = Eigen::Vector3d(1.,1.,1.);
    this->objects[0].noCheck = true;

    this->objects[1].filename = "./res/inner_box.obj";
    this->objects[1].min = Eigen::Vector3d(0.05, 0.3, 0.2);
    this->objects[1].max = Eigen::Vector3d(0.1, 0.5, 0.3);
    this->objects[1].noCheck = false;

    // simulation domain boundary
    this->sim_boundary_min = Eigen::Vector3d(-1,-1,-1);
    this->sim_boundary_max = Eigen::Vector3d(1, 1, 1);
    this->simbound_active = true;

    this->surface_reco_method = 1;  // 0: dense; 1: sparse
    this->pressure_solver_method = 0;   // 0: wcsph, 1: pbf
    this->n_iterations_pbf = 10;

    this->dt_default = 0.0005;
    this->t_between_frames = 0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.1;
    this->v_b = 0.0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
    this->assignment = "assignment5/boundary_wetting";

    this->emitters.resize(1);
    this->emitters[0].dir = {0,1,0};
    this->emitters[0].origin = {.06,0.15,.4};
    this->emitters[0].r = 0.1;
    this->emitters[0].velocity = 0.7;
    this->emitters[0].alternating = false;
    this->emitters[0].emission_freq = 1;
    this->emitters[0].emit_counter = 1000;    


    this->surface_tension = true;
    this->cohesion_coefficient = 0.0; //0.05;
    this->adhesion_coefficient = 0.1; //0.01;
}

void learnSPH::simulations_setup::Simulations::boundary_wetting_only_cohesion(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(0);
    this->fluid_end.resize(0);
    this->fluid_velocities.resize(0);

    this->objects.resize(2);
    this->objects[0].filename = "./res/boundary_cube.obj";
    this->objects[0].min = Eigen::Vector3d(-0.02,-0.02,-0.02);
    this->objects[0].max = Eigen::Vector3d(1.,1.,1.);
    this->objects[0].noCheck = true;

    this->objects[1].filename = "./res/inner_box.obj";
    this->objects[1].min = Eigen::Vector3d(0.05, 0.3, 0.2);
    this->objects[1].max = Eigen::Vector3d(0.1, 0.5, 0.3);
    this->objects[1].noCheck = false;

    // simulation domain boundary
    this->sim_boundary_min = Eigen::Vector3d(-1,-1,-1);
    this->sim_boundary_max = Eigen::Vector3d(1, 1, 1);
    this->simbound_active = true;

    this->surface_reco_method = 1;  // 0: dense; 1: sparse
    this->pressure_solver_method = 0;   // 0: wcsph, 1: pbf
    this->n_iterations_pbf = 10;

    this->dt_default = 0.0005;
    this->t_between_frames = 0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.1;
    this->v_b = 0.0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
    this->assignment = "assignment5/boundary_wetting";

    this->emitters.resize(1);
    this->emitters[0].dir = {0,1,0};
    this->emitters[0].origin = {.06,0.15,.4};
    this->emitters[0].r = 0.1;
    this->emitters[0].velocity = 0.7;
    this->emitters[0].alternating = false;
    this->emitters[0].emission_freq = 1;
    this->emitters[0].emit_counter = 1000;    


    this->surface_tension = true;
    this->cohesion_coefficient = 0.05; //0.05;
    this->adhesion_coefficient = 0.0; //0.01;
}

void learnSPH::simulations_setup::Simulations::boundary_wetting_cohesion_and_adhesion(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(0);
    this->fluid_end.resize(0);
    this->fluid_velocities.resize(0);

    this->objects.resize(2);
    this->objects[0].filename = "./res/boundary_cube.obj";
    this->objects[0].min = Eigen::Vector3d(-0.02,-0.02,-0.02);
    this->objects[0].max = Eigen::Vector3d(1.,1.,1.);
    this->objects[0].noCheck = true;

    this->objects[1].filename = "./res/inner_box.obj";
    this->objects[1].min = Eigen::Vector3d(0.05, 0.3, 0.2);
    this->objects[1].max = Eigen::Vector3d(0.1, 0.5, 0.3);
    this->objects[1].noCheck = false;

    // simulation domain boundary
    this->sim_boundary_min = Eigen::Vector3d(-1,-1,-1);
    this->sim_boundary_max = Eigen::Vector3d(1, 1, 1);
    this->simbound_active = true;

    this->surface_reco_method = 1;  // 0: dense; 1: sparse
    this->pressure_solver_method = 0;   // 0: wcsph, 1: pbf
    this->n_iterations_pbf = 10;

    this->dt_default = 0.0005;
    this->t_between_frames = 0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.0025;
    this->v_b = 0.0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
    this->assignment = "assignment5/boundary_wetting";

    this->emitters.resize(1);
    this->emitters[0].dir = {0,1,0};
    this->emitters[0].origin = {.1,0.15,.4};
    this->emitters[0].r = 0.1;
    this->emitters[0].velocity = 0.7;
    this->emitters[0].alternating = false;
    this->emitters[0].emission_freq = 1;
    this->emitters[0].emit_counter = 1000;    


    this->surface_tension = true;
    this->cohesion_coefficient = 0.05; //0.05;
    this->adhesion_coefficient = 0.1; //0.01;
}


void learnSPH::simulations_setup::Simulations::galton_board(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(0);
    this->fluid_velocities.resize(0);
    this->fluid_end.resize(0);

    // this->fluid_begin[0] = Eigen::Vector3d(0.054674, -0.18376, 1.5343);
    // this->fluid_end[0] = Eigen::Vector3d(-0.09, 0.18445, 1.7432);
    // this->fluid_velocities[0] = Eigen::Vector3d({0,0,0});

    this->objects.resize(3);

    this->objects[0].filename = "./res/galton-board-cylinders.obj";
    this->objects[1].filename = "./res/galton-board-front.obj";
    this->objects[2].filename = "./res/galton-board.obj";
    this->objects[0].noCheck = true;
    this->objects[1].noCheck = true;
    this->objects[2].noCheck = true;


    // simulation domain boundary
    this->sim_boundary_min = Eigen::Vector3d(-0.2,-0.5,0);
    this->sim_boundary_max = Eigen::Vector3d(0.2,0.5,2);
    this->simbound_active = true;

    this->surface_reco_method = 1;  // 0: dense; 1: sparse
    this->pressure_solver_method = 1;   // 0: wcsph, 1: pbf
    this->n_iterations_pbf = 10;

    this->dt_default = 0.001;
    this->t_between_frames = 0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.0025;
    this->v_b = 0.0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
    this->assignment = "final/galton_board";

    this->emitters.resize(2);
    this->emitters[0].dir = {0,-1,0};
    this->emitters[0].origin = {0.00157, 0.48307, 1.9};
    this->emitters[0].r = 0.035;
    this->emitters[0].velocity = 1;
    this->emitters[0].emission_freq = 1;
    this->emitters[0].emit_counter = 500;   

    this->emitters[1].dir = {0,1,0};
    this->emitters[1].origin = {0.00157, -0.48307, 1.9};
    this->emitters[1].r = 0.035;
    this->emitters[1].velocity = 1;
    this->emitters[1].emission_freq = 1;
    this->emitters[1].emit_counter = 500;   

    this->surface_tension = true;
    this->cohesion_coefficient = 0.05; //0.05;
    this->adhesion_coefficient = 0.01; //0.01;
}


void learnSPH::simulations_setup::Simulations::fountain_with_path(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(1);
    this->fluid_end.resize(1);
    this->fluid_velocities.resize(1);

    this->fluid_begin[0] = Eigen::Vector3d(-0.2, -0.2, -0.35);
    this->fluid_end[0] = Eigen::Vector3d(0.2, .2, -0.3);
    this->fluid_velocities[0] = Eigen::Vector3d(0.0, 0.0, 0.0);

    this->surface_reco_method = 1;  // 0: dense; 1: sparse
    this->pressure_solver_method = 1;   // 0: wcsph, 1: pbf
    this->n_iterations_pbf = 10;

    this->dt_default = 0.005;
    this->t_between_frames = 0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.0025;
    this->v_b = 0.0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.8);
    this->assignment = "assignment5/fountain_with_path";
    this->simTime = 5;

    this->surface_tension = true;
    this->cohesion_coefficient = 0.01;
    this->adhesion_coefficient = 0.01;

    this->emitters.resize(1);
    this->emitters[0].dir = {0,1,0};
    this->emitters[0].origin = {0.0,-0.2,-.2};
    this->emitters[0].r = 0.05;
    this->emitters[0].velocity = .7;
    this->emitters[0].alternating = false;
    this->emitters[0].emission_freq = 1;
    this->emitters[0].emit_counter = 1000;  

    this->objects.resize(1);
    this->objects[0].filename = "./res/boxes_connected.obj";
    this->objects[0].min = Eigen::Vector3d(-0.3, -0.3, -0.8);
    this->objects[0].max = Eigen::Vector3d(0.3, 2.0, 0.4);
    this->objects[0].noCheck = true;

    this->sim_boundary_min = this->objects[0].min;
    this->sim_boundary_max = this->objects[0].max;
    this->simbound_active = true;

}

void learnSPH::simulations_setup::Simulations::multiple_fountains_with_path(){
    this->particle_radius = 0.004;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(1);
    this->fluid_end.resize(1);
    this->fluid_velocities.resize(1);

    this->fluid_begin[0] = Eigen::Vector3d(-0.24, 0.02,-0.16);
    this->fluid_end[0] = Eigen::Vector3d(0.3, .1, -0.1);
    this->fluid_velocities[0] = Eigen::Vector3d(0.0, 0.0, 0.0);

    this->surface_reco_method = 1;  // 0: dense; 1: sparse
    this->pressure_solver_method = 1;   // 0: wcsph, 1: pbf
    this->n_iterations_pbf = 10;

    this->dt_default = 0.00025;
    this->t_between_frames = 0.004;
    this->B = 1000 * 1.02;
    this->v_f = 0.0025;
    this->v_b = 0.0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.8);
    this->assignment = "assignment5/maze";
    this->simTime = 20;

    this->surface_tension = true;
    this->cohesion_coefficient = 0.05;
    this->adhesion_coefficient = 0.01;

    this->emitters.resize(3);
    this->emitters[0].dir = {0,1,0};
    this->emitters[0].origin = {0.0,-0.1,0.};
    this->emitters[0].r = 0.024;
    this->emitters[0].velocity = 1.;
    this->emitters[0].alternating = false;
    this->emitters[0].emission_freq = 1;
    this->emitters[0].emit_counter = 10000;  

    this->emitters[1].dir = {0,1,0};
    this->emitters[1].origin = {0.15,-0.1,0.};
    this->emitters[1].r = 0.024;
    this->emitters[1].velocity = 1.;
    this->emitters[1].alternating = false;
    this->emitters[1].emission_freq = 1;
    this->emitters[1].emit_counter = 10000;  

    this->emitters[2].dir = {0,1,0};
    this->emitters[2].origin = {-0.15,-0.1,0.};
    this->emitters[2].r = 0.024;
    this->emitters[2].velocity = 1.;
    this->emitters[2].alternating = false;
    this->emitters[2].emission_freq = 1;
    this->emitters[2].emit_counter = 10000;  

    this->objects.resize(1);
    this->objects[0].filename = "./res/maze_small.obj";
    this->objects[0].min = Eigen::Vector3d(-.6, -0.15, -0.4);
    this->objects[0].max = Eigen::Vector3d(.6, 1.0, .1);
    this->objects[0].noCheck = true;

    this->sim_boundary_min = this->objects[0].min;
    this->sim_boundary_max = this->objects[0].max;
    this->simbound_active = true;
}


void learnSPH::simulations_setup::Simulations::simple_emitter_mid_cohesion(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(0);
    this->fluid_end.resize(0);
    this->fluid_velocities.resize(0);

    this->objects.resize(1);
    this->objects[0].filename = "./res/vessel.obj";
    this->objects[0].noCheck = true;

    // simulation domain boundary
    this->sim_boundary_min = Eigen::Vector3d(0,0,0);
    this->sim_boundary_max = Eigen::Vector3d(1, 1, .5);
    this->simbound_active = true;

    this->surface_reco_method = 1;  // 0: dense; 1: sparse
    this->pressure_solver_method = 1;   // 0: wcsph, 1: pbf
    this->n_iterations_pbf = 10;

    this->dt_default = 0.005;
    this->t_between_frames = 0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.0025;
    this->v_b = 0.0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
    this->assignment = "assignment5/simple_emitter_mid_cohesion";

    this->emitters.resize(1);
    this->emitters[0].dir = {0,1,0};
    this->emitters[0].origin = {.5,.1,.4};
    this->emitters[0].r = 0.05;
    this->emitters[0].velocity = 1;
    this->emitters[0].alternating = false;
    this->emitters[0].emission_freq = 1;
    this->emitters[0].emit_counter = 1000;

    this->surface_tension = true;
    this->cohesion_coefficient = 0.2;
    this->adhesion_coefficient = 0.0;
}



void learnSPH::simulations_setup::Simulations::simple_emitter_high_cohesion(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(0);
    this->fluid_end.resize(0);
    this->fluid_velocities.resize(0);

    this->objects.resize(1);
    this->objects[0].filename = "./res/vessel.obj";
    this->objects[0].noCheck = true;

    // simulation domain boundary
    this->sim_boundary_min = Eigen::Vector3d(0,0,0);
    this->sim_boundary_max = Eigen::Vector3d(1, 1, .5);
    this->simbound_active = true;

    this->surface_reco_method = 1;  // 0: dense; 1: sparse
    this->pressure_solver_method = 1;   // 0: wcsph, 1: pbf
    this->n_iterations_pbf = 10;

    this->dt_default = 0.005;
    this->t_between_frames = 0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.0025;
    this->v_b = 0.0;
    // this->v_f = 0.5;
    // this->v_b = 0.0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
    this->assignment = "assignment5/simple_emitter_high_cohesion";

    this->emitters.resize(1);
    this->emitters[0].dir = {0,1,0};
    this->emitters[0].origin = {.5,.1,.4};
    this->emitters[0].r = 0.05;
    this->emitters[0].velocity = 1;
    this->emitters[0].alternating = false;
    this->emitters[0].emission_freq = 1;
    this->emitters[0].emit_counter = 1000;

    this->surface_tension = true;
    this->cohesion_coefficient = 0.4;
    this->adhesion_coefficient = 0.0;
}



void learnSPH::simulations_setup::Simulations::water_droplet_no_adhesion(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(0);
    this->fluid_end.resize(0);
    this->fluid_velocities.resize(0);

    // this->fluid_begin[0] = Eigen::Vector3d(0.5, 0.2, 0.12);
    // this->fluid_end[0] = Eigen::Vector3d(0.6, 0.4, 0.22);
    // this->fluid_velocities[0] = Eigen::Vector3d(0.0, 0.0, 0.0);

    this->objects.resize(1);
    this->objects[0].filename = "./res/sphere.obj";
    this->objects[0].noCheck = true;

    // simulation domain boundary
    this->sim_boundary_min = Eigen::Vector3d(-0.5,-0.5,-0.5);
    this->sim_boundary_max = Eigen::Vector3d(1., 1., 1.);
    this->simbound_active = true;

    this->surface_reco_method = 1;  // 0: dense; 1: sparse
    this->pressure_solver_method = 1;   // 0: wcsph, 1: pbf
    this->n_iterations_pbf = 10;

    this->dt_default = 0.00025;
    this->t_between_frames = 0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.0025;
    this->v_b = 0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, 0.);
    this->assignment = "assignment5/water_droplet_no_adhesion";
    this->simTime = 5;

    this->surface_tension = false;
    this->cohesion_coefficient = 0.0;
    this->adhesion_coefficient = 0.0;

    this->emitters.resize(1);
    this->emitters[0].dir = {0,0,-1};
    this->emitters[0].origin = {0.0,0.0,.25};
    this->emitters[0].r = 0.025;
    this->emitters[0].velocity = 1;
    this->emitters[0].alternating = false;
    this->emitters[0].emission_freq = 1;
    this->emitters[0].emit_counter = 1000;
}


void learnSPH::simulations_setup::Simulations::water_droplet_mid_adhesion(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(0);
    this->fluid_end.resize(0);
    this->fluid_velocities.resize(0);

    // this->fluid_begin[0] = Eigen::Vector3d(0.5, 0.2, 0.12);
    // this->fluid_end[0] = Eigen::Vector3d(0.6, 0.4, 0.22);
    // this->fluid_velocities[0] = Eigen::Vector3d(0.0, 0.0, 0.0);

    this->objects.resize(1);
    this->objects[0].filename = "./res/sphere.obj";
    this->objects[0].noCheck = true;

    // simulation domain boundary
    this->sim_boundary_min = Eigen::Vector3d(-0.5,-0.5,-0.5);
    this->sim_boundary_max = Eigen::Vector3d(1., 1., 1.);
    this->simbound_active = true;

    this->surface_reco_method = 1;  // 0: dense; 1: sparse
    this->pressure_solver_method = 1;   // 0: wcsph, 1: pbf
    this->n_iterations_pbf = 10;

    this->dt_default = 0.00025;
    this->t_between_frames = 0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.0025;
    this->v_b = 0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, 0.);
    this->assignment = "assignment5/water_droplet_mid_adhesion";
    this->simTime = 5;

    this->surface_tension = true;
    this->cohesion_coefficient = 0.0;
    this->adhesion_coefficient = 1.0;

    this->emitters.resize(1);
    this->emitters[0].dir = {0,0,-1};
    this->emitters[0].origin = {0.0,0.0,.25};
    this->emitters[0].r = 0.025;
    this->emitters[0].velocity = 1;
    this->emitters[0].alternating = false;
    this->emitters[0].emission_freq = 1;
    this->emitters[0].emit_counter = 1000;
}


void learnSPH::simulations_setup::Simulations::water_droplet_high_adhesion(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(0);
    this->fluid_end.resize(0);
    this->fluid_velocities.resize(0);

    // this->fluid_begin[0] = Eigen::Vector3d(0.5, 0.2, 0.12);
    // this->fluid_end[0] = Eigen::Vector3d(0.6, 0.4, 0.22);
    // this->fluid_velocities[0] = Eigen::Vector3d(0.0, 0.0, 0.0);

    this->objects.resize(1);
    this->objects[0].filename = "./res/sphere.obj";
    this->objects[0].noCheck = true;

    // simulation domain boundary
    this->sim_boundary_min = Eigen::Vector3d(-0.5,-0.5,-0.5);
    this->sim_boundary_max = Eigen::Vector3d(1., 1., 1.);
    this->simbound_active = true;

    this->surface_reco_method = 1;  // 0: dense; 1: sparse
    this->pressure_solver_method = 1;   // 0: wcsph, 1: pbf
    this->n_iterations_pbf = 10;

    this->dt_default = 0.00025;
    this->t_between_frames = 0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.0025;
    this->v_b = 0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, 0.);
    this->assignment = "assignment5/water_droplet_high_adhesion";
    this->simTime = 5;

    this->surface_tension = true;
    this->cohesion_coefficient = 0.0;
    this->adhesion_coefficient = 10.0;

    this->emitters.resize(1);
    this->emitters[0].dir = {0,0,-1};
    this->emitters[0].origin = {0.0,0.0,.25};
    this->emitters[0].r = 0.025;
    this->emitters[0].velocity = 1;
    this->emitters[0].alternating = false;
    this->emitters[0].emission_freq = 1;
    this->emitters[0].emit_counter = 1000;
}


void learnSPH::simulations_setup::Simulations::water_droplet_on_surface(){

}


void learnSPH::simulations_setup::Simulations::droplets_on_leaf(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->objects.resize(2);
    this->objects[0].filename = "./res/leaf_big.obj";
    this->objects[0].min = Eigen::Vector3d(-1., -1., -1.);
    this->objects[0].max = Eigen::Vector3d(1., 1., 1.2);
    this->objects[0].noCheck = true;

    this->objects[1].filename = "./res/floor2.obj";
    this->objects[1].min = Eigen::Vector3d(-2., -2., -1.1);
    this->objects[1].max = Eigen::Vector3d(2., 2., -0.9);
    this->objects[1].noCheck = true;

    this->sim_boundary_min = this->objects[0].min;
    this->sim_boundary_max = this->objects[0].max;
    this->simbound_active = true;

    this->surface_tension = true;
    this->cohesion_coefficient = 0.01;
    this->adhesion_coefficient = 1.;

    this->surface_reco_method = 1;  // 0: dense; 1: sparse
    this->pressure_solver_method = 1;   // 0: wcsph, 1: pbf
    this->n_iterations_pbf = 5;

    this->dt_default = 0.0002;
    this->t_between_frames = 0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.1;
    this->v_b = 0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.8);
    this->assignment = "assignment5/droplets_on_leaf";
    this->simTime = 5;

    this->emitters.resize(20);
    this->emitters[0].dir = {0,0,-1};
    this->emitters[0].origin = {0.23,-0.03,0.4};
    this->emitters[0].r = 0.01;
    this->emitters[0].velocity = 1.;
    this->emitters[0].alternating = false;
    this->emitters[0].emission_freq = 4;
    this->emitters[0].emit_counter = 1000;  

    this->emitters[1].dir = {0,0,-1};
    this->emitters[1].origin = {0.15,-0.1,0.41};
    this->emitters[1].r = 0.01;
    this->emitters[1].velocity = 1.;
    this->emitters[1].alternating = false;
    this->emitters[1].emission_freq = 4;
    this->emitters[1].emit_counter = 1000;  

    this->emitters[2].dir = {0,0,-1};
    this->emitters[2].origin = {0.3,0.01,0.43};
    this->emitters[2].r = 0.01;
    this->emitters[2].velocity = 1.;
    this->emitters[2].alternating = false;
    this->emitters[2].emission_freq = 4;
    this->emitters[2].emit_counter = 1000;  

    this->emitters[3].dir = {0,0,-1};
    this->emitters[3].origin = {0.2,-0.01,0.42};
    this->emitters[3].r = 0.01;
    this->emitters[3].velocity = 1.;
    this->emitters[3].alternating = false;
    this->emitters[3].emission_freq = 4;
    this->emitters[3].emit_counter = 1000;  

    this->emitters[4].dir = {0,0,-1};
    this->emitters[4].origin = {0.1,0.005,0.45};
    this->emitters[4].r = 0.01;
    this->emitters[4].velocity = 1.;
    this->emitters[4].alternating = false;
    this->emitters[4].emission_freq = 3;
    this->emitters[4].emit_counter = 1000;  

    this->emitters[5].dir = {0,0,-1};
    this->emitters[5].origin = {0.0,0.002,0.4};
    this->emitters[5].r = 0.01;
    this->emitters[5].velocity = 1.;
    this->emitters[5].alternating = false;
    this->emitters[5].emission_freq = 4;
    this->emitters[5].emit_counter = 1000;  

    this->emitters[6].dir = {0,0,-1};
    this->emitters[6].origin = {0.05,0.1,0.4};
    this->emitters[6].r = 0.01;
    this->emitters[6].velocity = 1.;
    this->emitters[6].alternating = false;
    this->emitters[6].emission_freq = 3;
    this->emitters[6].emit_counter = 1000;  

    this->emitters[7].dir = {0,0,-1};
    this->emitters[7].origin = {0.13,0.11,0.42};
    this->emitters[7].r = 0.01;
    this->emitters[7].velocity = 1.;
    this->emitters[7].alternating = false;
    this->emitters[7].emission_freq = 4;
    this->emitters[7].emit_counter = 1000;  

    this->emitters[8].dir = {0,0,-1};
    this->emitters[8].origin = {0.2,0.09,0.44};
    this->emitters[8].r = 0.01;
    this->emitters[8].velocity = 1.;
    this->emitters[8].alternating = false;
    this->emitters[8].emission_freq = 5;
    this->emitters[8].emit_counter = 1000;  

    this->emitters[9].dir = {0,0,-1};
    this->emitters[9].origin = {0.27,0.08,0.43};
    this->emitters[9].r = 0.01;
    this->emitters[9].velocity = 1.;
    this->emitters[9].alternating = false;
    this->emitters[9].emission_freq = 4;
    this->emitters[9].emit_counter = 1000;  

    this->emitters[10].dir = {0,0,-1};
    this->emitters[10].origin = {0.04,-0.13,0.43};
    this->emitters[10].r = 0.01;
    this->emitters[10].velocity = 1.;
    this->emitters[10].alternating = false;
    this->emitters[10].emission_freq = 4;
    this->emitters[10].emit_counter = 1000; 

    this->emitters[11].dir = {0,0,-1};
    this->emitters[11].origin = {-0.1,-0.01,0.42};
    this->emitters[11].r = 0.01;
    this->emitters[11].velocity = 1.;
    this->emitters[11].alternating = false;
    this->emitters[11].emission_freq = 3;
    this->emitters[11].emit_counter = 1000;  

    this->emitters[12].dir = {0,0,-1};
    this->emitters[12].origin = {-0.03,0.07,0.4};
    this->emitters[12].r = 0.01;
    this->emitters[12].velocity = 1.;
    this->emitters[12].alternating = false;
    this->emitters[12].emission_freq = 5;
    this->emitters[12].emit_counter = 1000;   



    this->emitters[13].dir = {0,0,-1};
    this->emitters[13].origin = {0.01,-0.04,0.32};
    this->emitters[13].r = 0.01;
    this->emitters[13].velocity = 1.;
    this->emitters[13].alternating = false;
    this->emitters[13].emission_freq = 4;
    this->emitters[13].emit_counter = 1000;  

    this->emitters[14].dir = {0,0,-1};
    this->emitters[14].origin = {-0.1,-0.035,0.32};
    this->emitters[14].r = 0.01;
    this->emitters[14].velocity = 1.;
    this->emitters[14].alternating = false;
    this->emitters[14].emission_freq = 4;
    this->emitters[14].emit_counter = 1000;  

    this->emitters[15].dir = {0,0,-1};
    this->emitters[15].origin = {-0.3,-0.028,0.3};
    this->emitters[15].r = 0.01;
    this->emitters[15].velocity = 1.;
    this->emitters[15].alternating = false;
    this->emitters[15].emission_freq = 4;
    this->emitters[15].emit_counter = 1000;  

    this->emitters[16].dir = {0,0,-1};
    this->emitters[16].origin = {-0.1,-0.1,0.42};
    this->emitters[16].r = 0.01;
    this->emitters[16].velocity = 1.;
    this->emitters[16].alternating = false;
    this->emitters[16].emission_freq = 4;
    this->emitters[16].emit_counter = 1000;  

    this->emitters[17].dir = {0,0,-1};
    this->emitters[17].origin = {-0.05,-0.03,0.45};
    this->emitters[17].r = 0.01;
    this->emitters[17].velocity = 1.;
    this->emitters[17].alternating = false;
    this->emitters[17].emission_freq = 3;
    this->emitters[17].emit_counter = 1000;  

    this->emitters[18].dir = {0,0,-1};
    this->emitters[18].origin = {-0.26,-0.042,0.4};
    this->emitters[18].r = 0.01;
    this->emitters[18].velocity = 1.;
    this->emitters[18].alternating = false;
    this->emitters[18].emission_freq = 4;
    this->emitters[18].emit_counter = 1000;  

    this->emitters[19].dir = {0,0,-1};
    this->emitters[19].origin = {-0.05,-0.1,0.4};
    this->emitters[19].r = 0.01;
    this->emitters[19].velocity = 1.;
    this->emitters[19].alternating = false;
    this->emitters[19].emission_freq = 3;
    this->emitters[19].emit_counter = 1000;  

}



void learnSPH::simulations_setup::Simulations::single_flow(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(1);
    this->fluid_end.resize(1);
    this->fluid_velocities.resize(1);

    this->fluid_begin[0] = Eigen::Vector3d(0., 0., 0.);
    this->fluid_end[0] = Eigen::Vector3d(0.15, 0.78, 0.1);
    this->fluid_velocities[0] = Eigen::Vector3d(0.0, 0.0, 0.0);

    this->objects.resize(1);
    this->objects[0].filename = "./res/boundary.obj";
    this->objects[0].min = Eigen::Vector3d(-0.02, -0.02, -0.02);
    this->objects[0].max = Eigen::Vector3d(0.17, 0.8, 1.0);
    this->objects[0].noCheck = true;

    this->sim_boundary_min = Eigen::Vector3d(-1., -1., -1.);
    this->sim_boundary_max = Eigen::Vector3d(1., 1., 1.);
    this->simbound_active = true;

    this->surface_tension = false;
    this->cohesion_coefficient = 0.01;
    this->adhesion_coefficient = 1.;

    this->surface_reco_method = 1;  // 0: dense; 1: sparse
    this->pressure_solver_method = 1;   // 0: wcsph, 1: pbf
    this->n_iterations_pbf = 5;

    this->dt_default = 0.0002;
    this->t_between_frames = 0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.1;
    this->v_b = 0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.8);
    this->assignment = "assignment5/single_flow";
    this->simTime = 5;

    this->emitters.resize(1);
    this->emitters[0].dir = {0,0,-1};
    this->emitters[0].origin = {0.1,0.2,0.3};
    this->emitters[0].r = 0.01;
    this->emitters[0].velocity = 1.;
    this->emitters[0].alternating = false;
    this->emitters[0].emission_freq = 3;
    this->emitters[0].emit_counter = 1000;    

}


void learnSPH::simulations_setup::Simulations::puzzle(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->objects.resize(1);
    this->objects[0].filename = "./res/puzzle.obj";
    this->objects[0].min = Eigen::Vector3d(-0.4, -1., -1.);
    this->objects[0].max = Eigen::Vector3d(0.4, 1., 1.0);
    this->objects[0].noCheck = true;

    this->sim_boundary_min = this->objects[0].min;
    this->sim_boundary_max = this->objects[0].max;
    this->simbound_active = true;

    this->surface_tension = false;
    this->cohesion_coefficient = 0.01;
    this->adhesion_coefficient = 1.;

    this->surface_reco_method = 1;  // 0: dense; 1: sparse
    this->pressure_solver_method = 1;   // 0: wcsph, 1: pbf
    this->n_iterations_pbf = 5;

    this->dt_default = 0.0002;
    this->t_between_frames = 0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.1;
    this->v_b = 0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.8);
    this->assignment = "assignment5/puzzle";
    this->simTime = 10;

    this->emitters.resize(2);
    this->emitters[0].dir = {0,0,-1};
    this->emitters[0].origin = {0.1,-0.27,0.6};
    this->emitters[0].r = 0.025;
    this->emitters[0].velocity = 1.;
    this->emitters[0].alternating = false;
    this->emitters[0].emission_freq = 1;
    this->emitters[0].emit_counter = 2000; 

    this->emitters[1].dir = {0,0,-1};
    this->emitters[1].origin = {0.1,0.21,0.6};
    this->emitters[1].r = 0.025;
    this->emitters[1].velocity = 1.;
    this->emitters[1].alternating = false;
    this->emitters[1].emission_freq = 1;
    this->emitters[1].emit_counter = 2000; 
}