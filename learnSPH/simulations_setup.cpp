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

    this->sim_boundary_min = this->objects[0].min;
    this->sim_boundary_max = this->objects[0].max;
    this->simbound_active = true;

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

    this->dt_default = 0.001;
    this->t_between_frames =  0.008;
    this->B = 1000 * 1.02;
    this->v_f = 0.0025;
    this->v_b = 0.0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
    this->assignment = "assignment4/dam_break";
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
    this->n_iterations_pbf = 5;

    this->dt_default = 0.00025;
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
    this->emitters[0].origin = {.5,.1,.2};
    this->emitters[0].r = 0.1;
    this->emitters[0].velocity = .8;
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
    this->emitters[0].r = 0.1;
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
    this->pressure_solver_method = 0;   // 0: wcsph, 1: pbf

    this->dt_default = 0.0001;
    this->t_between_frames = 0.1;
    this->B = 1000 * 1.02;
    this->v_f = 0.1;
    this->v_b = 0;
    this->gravity = Eigen::Vector3d(0.0, 0.0, 0.0);
    this->assignment = "assignment5/water_droplet";

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
    this->adhesion_coefficient = 0.1; //0.01;
}


void learnSPH::simulations_setup::Simulations::galton_board(){
    this->particle_radius = 0.005;
    this->fluid_rest_density = 1000.0;

    this->fluid_begin.resize(0);
    this->fluid_end.resize(0);
    this->fluid_velocities.resize(0);

    this->objects.resize(16);
    this->objects[0].filename = "./res/galton_boundary.obj";
    this->objects[0].min = Eigen::Vector3d(0.0,0.0,0.0);
    this->objects[0].max = Eigen::Vector3d(0.5,1.1,3.0);
    this->objects[0].noCheck = true;

    this->objects[1].filename = "./res/galton_wall1.obj";
    this->objects[1].noCheck = false;
    this->objects[2].filename = "./res/galton_wall2.obj";
    this->objects[2].noCheck = false;
    this->objects[3].filename = "./res/galton_wall3.obj";
    this->objects[3].noCheck = false;
    this->objects[4].filename = "./res/galton_wall4.obj";
    this->objects[4].noCheck = false;

    this->objects[5].filename = "./res/galton_middle_pin1.obj";
    this->objects[5].noCheck = false;
    this->objects[6].filename = "./res/galton_middle_pin2.obj";
    this->objects[6].noCheck = false;
    this->objects[7].filename = "./res/galton_middle_pin3.obj";
    this->objects[7].noCheck = false;

    this->objects[8].filename = "./res/galton_middle+1_pin1.obj";
    this->objects[8].noCheck = false;
    this->objects[9].filename = "./res/galton_middle+1_pin2.obj";
    this->objects[9].noCheck = false;
    this->objects[10].filename = "./res/galton_middle+1_pin3.obj";
    this->objects[10].noCheck = false;

    this->objects[11].filename = "./res/galton_middle+2_pin1.obj";
    this->objects[11].noCheck = false;
    this->objects[12].filename = "./res/galton_middle+2_pin2.obj";
    this->objects[12].noCheck = false;

    this->objects[13].filename = "./res/galton_middle+3_pin1.obj";
    this->objects[13].noCheck = false;
    this->objects[14].filename = "./res/galton_middle+3_pin2.obj";
    this->objects[14].noCheck = false;

    this->objects[15].filename = "./res/galton_middle+4_pin1.obj";
    this->objects[15].noCheck = false;


    // simulation domain boundary
    this->sim_boundary_min = Eigen::Vector3d(-0.02,-0.02,-0.02);
    this->sim_boundary_max = Eigen::Vector3d(3.12, 3.12, 3.12);
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
    this->assignment = "assignment5/galton_board";

    this->emitters.resize(1);
    this->emitters[0].dir = {0,0,-1};
    this->emitters[0].origin = {.25,0.55,1.1};
    this->emitters[0].r = 0.1;
    this->emitters[0].velocity = 0.7;
    this->emitters[0].alternating = false;
    this->emitters[0].emission_freq = 1;
    this->emitters[0].emit_counter = 1000;   

    this->surface_tension = false;
    this->cohesion_coefficient = 0.0; //0.05;
    this->adhesion_coefficient = 0.0; //0.01;
}
