#include <array>
#include <cstdlib>
#include <iostream>


#include "marching_cubes.h"
#include "marching_cubes_lut.h"
#include "utils.h"

learnSPH::surface::MarchingCubes::MarchingCubes(double cellWidth, uint n_x, uint n_y, uint n_z,
                                                Eigen::Vector3d origin,
                                                learnSPH::types::ImplicitSurface implicitSurfaceFunction)
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

    this->levelSet.resize(n_vx * n_vy * n_vz);
    this->gridVertices.rehash(n_vx * n_vy * n_vz);
    this->triangles.rehash(n_ex * n_ey * n_ez);

    this->implicitSurfaceFunction = implicitSurfaceFunction;
}

void learnSPH::surface::MarchingCubes::get_Isosurface()
{
    uint vertexIdx, edgeIdx, edgeZIndex, intersecId, triangleId;
    Eigen::Vector3d offset_cords;
    std::array<bool, 8> levelSet;
    std::array<double, 8> vertex_signs;
    std::array<std::array<int, 3>, 5> triangulation;
    std::array<int, 2> vertexPair;
    std::array<Eigen::Vector3d, 8> vertices_coords;
    double intersecPoint, alpha, offset;

    intersecId = 0;
    triangleId = 0;

    // complete loop for finding the intersection points
    // Very high possibility for using wrong order of indices
    for (int i = 0; i < this->n_vx; ++i) {
        for (int j = 0; j < this->n_vy; ++j) {
            for (int k = 0; k < this->n_vz; ++k) {
                // saving vertex id to find edgeIDX for hasmap edge -> intersection
                vertexIdx = i * this->n_vy * this->n_vz + j * this->n_vz + k;
                // get level set and
                for (int l = 0; l < 8; ++l) {
                    // get real grid coordinates of vertices
                    vertices_coords[l] = {(CELL_VERTICES[l][0] + i) * cellWidth + origin.x(),
                                          (CELL_VERTICES[l][1] + j) * cellWidth + origin.y(),
                                          (CELL_VERTICES[l][2] + k) * cellWidth + origin.z()};
                    // compute sign distance value with implicit function
                    vertex_signs[l] = this->implicitSurfaceFunction(vertices_coords[l]);
                    // get levelset by simply checking if vertex inside implicit surface
                    levelSet[l] = vertex_signs[l] < 0;
                }
                // get corresponding triangulation for current cell
                triangulation = get_marching_cubes_cell_triangulation(levelSet);
                // assign first triple of triangle and immediatly skip the rest if first value is -1
                std::array<int, 3> triangle = triangulation[0];
                while (triangle[0] >= 0) {
                    // each entry in triangle corresponds to an edge in the cell
                    for (int l = 0; l < 3; ++l) {
                        // get vertex pair for first edge
                        vertexPair = CELL_EDGES[triangle[l]];
                        // compute absolute edge id
                        edgeIdx = vertexIdx * 3 + CELL_EDGES_DIRECTION[triangle[l]];
                        // compute intersection point
                        alpha = vertex_signs[vertexPair[0]] /
                                (vertex_signs[vertexPair[0]] - vertex_signs[vertexPair[1]]);
                        offset = (1.0 - alpha) * vertex_signs[vertexPair[0]] +
                                 alpha * vertex_signs[vertexPair[1]];
                        // intersection point is added to first vertex of vertices pair to get the
                        // real space vertex for intersection point
                        offset_cords = vertices_coords[vertexPair[0]];
                        offset_cords[CELL_EDGES_DIRECTION[triangle[l]]] += offset;
                        // put in the vector
                        this->intersections[intersecId] = offset_cords;
                        // add (key, value) pair to hashmap
                        this->edgeIntersection[edgeIdx] = intersecId;
                        this->triangles[triangleId][l]  = edgeIdx;
                        ++intersecId;
                    }
                    ++triangleId;
                }
            }
        }
    }
}

void learnSPH::surface::MarchingCubes::compute_normals(){
    if (this->triangles.size() == 0){
        std::cout << "No polygons to compute normals for!" << "\n";
        exit(-1);
    }
    

}

