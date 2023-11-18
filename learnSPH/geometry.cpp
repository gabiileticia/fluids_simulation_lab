#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <vector>

#include "../learnSPH/io.h"
#include "../learnSPH/sampling.h"
#include "../learnSPH/geometry.h"


void learnSPH::geometry::load_n_sample_boundary(std::vector<Eigen::Vector3d>& output
                            , std::vector<std::string> file_name
                            , double boundary_sampling_distance){
	std::vector<Eigen::Vector3d> aux_positions_vector;
	for (int i=0; i<file_name.size(); ++i){
		const std::vector<learnSPH::TriMesh> boundary_meshes = learnSPH::read_tri_meshes_from_obj(file_name[i]);
		const learnSPH::TriMesh& boundary = boundary_meshes[0];

		learnSPH::sampling::triangle_mesh(output
										, boundary.vertices
										, boundary.triangles
										, boundary_sampling_distance);
	}
}

void learnSPH::geometry::load_n_sample_fluids(std::vector<Eigen::Vector3d>& output_positions
                            , std::vector<Eigen::Vector3d>& output_velocities
                            , std::vector<Eigen::Vector3d> fluid_begin
                            , std::vector<Eigen::Vector3d> fluid_end
                            , double fluid_sampling_distance
							, std::vector<Eigen::Vector3d> fluid_velocities)
{
	std::vector<Eigen::Vector3d> aux_positions_vector;
	std::vector<Eigen::Vector3d> aux_velocities_vector;
	int begin_vector = 0;
	for (int i=0; i<fluid_begin.size(); ++i){

		learnSPH::sampling::fluid_box(output_positions, fluid_begin[i], fluid_end[i],fluid_sampling_distance);

		for(int j = begin_vector; j < output_positions.size();++j){
			output_velocities.push_back(fluid_velocities[i]);
		}
		begin_vector = output_velocities.size();
	}
}
    

