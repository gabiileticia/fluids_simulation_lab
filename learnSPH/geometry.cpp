#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <Eigen/Dense>

#include "../learnSPH/io.h"
#include "../learnSPH/sampling.h"
#include "../learnSPH/geometry.h"


void learnSPH::geometry::load_n_sample_boundary(std::vector<Eigen::Vector3d>& output
                            , std::string file_name
                            , double boundary_sampling_distance){
    const std::vector<learnSPH::TriMesh> boundary_meshes = learnSPH::read_tri_meshes_from_obj(file_name);
	const learnSPH::TriMesh& boundary = boundary_meshes[0];

	learnSPH::sampling::triangle_mesh(output
									, boundary.vertices
									, boundary.triangles
									, boundary_sampling_distance);
}