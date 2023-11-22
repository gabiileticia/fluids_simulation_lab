#include <array>
#include <cstdlib>
#include <iostream>

#include "marching_cubes.h"
#include "marching_cubes_lut.h"
#include "utils.h"

learnSPH::surface::MarchingCubes::MarchingCubes(
    double cellWidth, uint n_x, uint n_y, uint n_z, Eigen::Vector3d origin,
    learnSPH::types::ImplicitSurface implicitSurfaceFunction, void *funcArgs, double epsilon)
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

    this->levelSet.resize(n_vx * n_vy * n_vz);
    this->triangles.resize(n_cx * n_cy * n_cz * 5);
    this->intersections.resize(this->triangles.size() * 3);

    this->implicitSurfaceFunction = implicitSurfaceFunction;
    this->funcArgs                = funcArgs;
}

void learnSPH::surface::MarchingCubes::get_Isosurface()
{
    uint vertexIdx, edgeIdx, edgeZIndex, intersecId, triangleId;
    Eigen::Vector3d intersection;
    std::array<bool, 8> levelSet;
    std::array<double, 8> vertex_signs;
    std::array<std::array<int, 3>, 5> triangulation;
    std::array<int, 2> vertexPair;
    std::array<Eigen::Vector3d, 8> vertices_coords;
    double intersecPoint, alpha;

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
                    vertex_signs[l] =
                        this->implicitSurfaceFunction(vertices_coords[l], this->funcArgs);
                    // get levelset by simply checking if vertex inside implicit surface
                    levelSet[l] = vertex_signs[l] > 0;
                    // if (levelSet[l] == false) {
                    //     std::cout << "Found something!"
                    //               << "\n";
                    // }
                }
                // get corresponding triangulation for current cell
                triangulation = get_marching_cubes_cell_triangulation(levelSet);
                // assign first triple of triangle and immediatly skip the rest if first value is -1
                for (int m = 0; m < 5; ++m) {
                    std::array<int, 3> triangle = triangulation[m];
                    if (triangle[0] >= 0) {
                        // each entry in triangle corresponds to an edge in the cell
                        for (int l = 0; l < 3; ++l) {
                            // get vertex pair for first edge
                            vertexPair = CELL_EDGES[triangle[l]];
                            // compute absolute edge id
                            vertexIdx = learnSPH::utils::cubeVertex2VertexIndex(
                                vertexIdx, vertexPair[0], this->n_vx, this->n_vy, this->n_vz);
                            edgeIdx = vertexIdx * 3 + CELL_EDGES_DIRECTION[triangle[l]];
                            // compute intersection point
                            alpha = vertex_signs[vertexPair[0]] /
                                    (vertex_signs[vertexPair[0]] - vertex_signs[vertexPair[1]]);
                            intersection = (1.0 - alpha) * vertices_coords[vertexPair[0]] +
                                     alpha * vertices_coords[vertexPair[1]];
                            // put in the vector
                            this->intersections[intersecId] = intersection;
                            // add (key, value) pair to hashmap
                            this->edgeIntersection[edgeIdx] = intersecId;
                            this->triangles[triangleId][l]  = edgeIdx;
                            ++intersecId;
                        }
                        ++triangleId;
                    } else {
                        break;
                    }
                }
            }
        }
    }
    std::cout << "Finished retrieving ISOSurface!";
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
