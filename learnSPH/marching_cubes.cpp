#include <array>
#include <cstdlib>
#include <iostream>

#include "marching_cubes.h"
#include "marching_cubes_lut.h"
#include "utils.h"


learnSPH::surface::MarchingCubes::MarchingCubes(
    double cellWidth, uint n_x, uint n_y, uint n_z, Eigen::Vector3d origin, double epsilon, bool implicitFlag)
{
    this->origin = origin;
    this->n_cx   = n_x;
    this->n_cy   = n_y;
    this->n_cz   = n_z;
    this->n_vx   = n_x + 1;
    this->n_vy   = n_y + 1;
    this->n_vz   = n_z + 1;
    this->n_ex   = n_x;
    this->n_ey   = n_y;
    this->n_ez   = n_z;
    // Currently requires that realspace grid axis are divisible by the cellWidth \
            May lead to problems later.
    this->cellWidth = cellWidth;
    this->epsilon   = epsilon;

    this->triangles.resize(n_cx * n_cy * n_cz * 5);
    this->intersections.resize(this->triangles.size() * 3);
    this->implicitFlag = implicitFlag;
}

void learnSPH::surface::MarchingCubes::compute_normals()
{
    int a, b, c;
    Eigen::Vector3d v_a, v_b, v_c;
    if (this->triangles.size() == 0) {
        std::cout << "No polygons to compute normals for!"
                  << "\n";
        exit(-1);
    }
    this->intersectionNormals.resize(this->intersections.size());

    for (int i = 0; i < this->triangles.size(); ++i) {
        // get vertex indices
        a = this->edgeIntersection[triangles[i][0]];
        b = this->edgeIntersection[triangles[i][1]];
        c = this->edgeIntersection[triangles[i][2]];
        // get vertex values
        v_a = this->intersections[a];
        v_b = this->intersections[b];
        v_c = this->intersections[c];

        // compute normals with finite difference
        using namespace learnSPH::utils;
        this->intersectionNormals[a] =
            implicitVertexNormal(this->implicitSurfaceFunction, v_a, this->epsilon, this->funcArgs);
        this->intersectionNormals[b] =
            implicitVertexNormal(this->implicitSurfaceFunction, v_b, this->epsilon, this->funcArgs);
        this->intersectionNormals[c] =
            implicitVertexNormal(this->implicitSurfaceFunction, v_c, this->epsilon, this->funcArgs);
    }
}

void learnSPH::surface::MarchingCubes::compute_normals_gl(std::vector<double> &level_set)
{

    if (this->intersections.size() == 0) {
        std::cout << "No mesh vertices to compute normals for!"
                  << "\n";
        exit(-1);
    }
    // Unit vectors
    std::vector<Eigen::Vector3d> directions(4);
    directions[0] = Eigen::Vector3d(1, 0, 0);
    directions[1] = Eigen::Vector3d(0, 1, 0);
    directions[2] = Eigen::Vector3d(0, 0, 1);
    directions[3] = Eigen::Vector3d(1, 0, 0);

    Eigen::Vector3d fin_diff;
    Eigen::Vector3d normal;

    static Eigen::Vector3d ex = Eigen::Vector3d(1, 0, 0);
    static Eigen::Vector3d ey = Eigen::Vector3d(0, 1, 0);
    static Eigen::Vector3d ez = Eigen::Vector3d(0, 0, 1);

    this->intersectionNormals.resize(this->intersections.size(), Eigen::Vector3d(0, 0, 0));
    uint origin_vertex, end_vertex;
    int direction;
    std::vector<int> sum_triangles(intersectionNormals.size(), 0);

    for(int i=0; i< triangles.size(); i++){
        for (int j=0;j<3;j++){

            Eigen::Vector3d u = intersections[triangles[i][j]] - intersections[triangles[i][(j + 1)%3]];
            Eigen::Vector3d v = intersections[triangles[i][j]] - intersections[triangles[i][(j + 2)%3]];

            normal[0] = u.y() * v.z() - u.z() * v.y(); 
            normal[1] = u.z() * v.x() - u.x() * v.z();
            normal[2] = u.x() * v.y() - u.y() * v.x();

            this->intersectionNormals[triangles[i][j]] += normal.normalized();
            sum_triangles[triangles[i][j]] += 1; 

        }
    }
    for(int i=0; i<intersectionNormals.size();i++){
        intersectionNormals[i] = (intersectionNormals[i] / sum_triangles[i]).normalized();
    }
}


void learnSPH::surface::MarchingCubes::get_isosurface(std::vector<double> &level_set)
{
    uint vertex_id, vertex_id_end, edge_id, edgeZIndex, intersecId, triangleId;
    std::array<uint, 8> vertex_ids;
    Eigen::Vector3d intersection;
    std::array<bool, 8> vertex_signs;
    std::array<std::array<int, 3>, 5> triangulation;
    std::array<int, 2> vertex_pair;
    std::array<Eigen::Vector3d, 8> vertices_coords;
    double intersecPoint, alpha;

    intersecId = 0;
    triangleId = 0;

    std::vector<Eigen::Vector3d> directions(3);
    directions[0] = Eigen::Vector3d(1, 0, 0);
    directions[1] = Eigen::Vector3d(0, 1, 0);
    directions[2] = Eigen::Vector3d(0, 0, 1);

    for(int i=0; i < this->n_cx; i++){
        for(int j=0; j < this->n_cy; j++){
            for(int k=0; k < this->n_cz; k++){
        
                for (int l = 0; l < 8; l++) {
                    vertex_ids[l] = (i + CELL_VERTICES[l][0]) * this->n_vy * this->n_vz + (j + CELL_VERTICES[l][1]) * this->n_vz + (k + CELL_VERTICES[l][2]);
                    vertex_signs[l] = level_set[vertex_ids[l]] < 0;   
                    // get real grid coordinates of vertices
                    vertices_coords[l] = {(CELL_VERTICES[l][0] + i) * cellWidth + origin.x(),
                                          (CELL_VERTICES[l][1] + j) * cellWidth + origin.y(),
                                          (CELL_VERTICES[l][2] + k) * cellWidth + origin.z()}; 
                }
                triangulation = get_marching_cubes_cell_triangulation(vertex_signs);
                for(int m=0 ;m<5 ; m++){
                    std::array<int, 3> triangle = triangulation[m];
                    if (triangle[0] >= 0) {
                        for (int l = 0; l < 3; ++l) {
                            vertex_pair = CELL_EDGES[triangle[l]];
                            edge_id = vertex_ids[vertex_pair[0]] * 3 + CELL_EDGES_DIRECTION[triangle[l]];
                            alpha = level_set[vertex_ids[vertex_pair[0]]] /
                                    (level_set[vertex_ids[vertex_pair[0]]] - level_set[vertex_ids[vertex_pair[1]]]);
                            intersection = (1.0 - alpha) * vertices_coords[vertex_pair[0]] +
                                           alpha * vertices_coords[vertex_pair[1]];
                            this->intersections[intersecId] = intersection;
                            this->edgeIntersection[intersecId] = edge_id;
                            this->triangles[triangleId][l] = intersecId;
                            intersecId++;
                        }
                        triangleId++;
                    } else {
                        break;
                    }
                }
            }
        }
    }
    this->triangles.resize(triangleId*3);
    this->intersections.resize(intersecId);
}
