#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sys/types.h>
#include <system_error>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "marching_cubes.h"
#include "marching_cubes_lut.h"
#include "utils.h"

learnSPH::surface::MarchingCubes::MarchingCubes(
    double cellWidth, uint n_x, uint n_y, uint n_z, Eigen::Vector3d origin,
    learnSPH::theta_functions::FluidThetaFunction &thetaFunction, double epsilon, bool implicitFlag,
    double c)
    : origin(origin), n_cx(n_x), n_cy(n_y), n_cz(n_z), n_vx(n_x + 1), n_vy(n_y + 1), n_vz(n_z + 1),
      n_ex(n_x), n_ey(n_y), n_ez(n_z), cellWidth(cellWidth), epsilon(epsilon),
      thetaFunction(thetaFunction), c(c)
{
    this->triangles.resize(n_cx * n_cy * n_cz * 5);
    this->intersections.resize(this->triangles.size() * 3);

    this->implicitFlag = implicitFlag;
}

void learnSPH::surface::MarchingCubes::get_Isosurface(std::vector<double> &level_set)
{
    uint vertexIdx, edgeIdx, edgeZIndex, intersecId, triangleId;
    Eigen::Vector3d intersection;
    std::array<bool, 8> vertex_signs;
    std::array<uint, 8> vertex_ids;
    std::array<std::array<int, 3>, 5> triangulation;
    std::array<int, 2> vertexPair;
    std::array<Eigen::Vector3d, 8> vertices_coords;
    double intersecPoint, alpha;

    intersecId = 0;
    triangleId = 0;

    // complete loop for finding the intersection points
    // Very high possibility for using wrong order of indices
    for (int i = 0; i < this->n_cx; ++i) {
        for (int j = 0; j < this->n_cy; ++j) {
            for (int k = 0; k < this->n_cz; ++k) {
                for (int l = 0; l < 8; ++l) {
                    // get real grid coordinates of vertices
                    vertex_ids[l] = (i + CELL_VERTICES[l][0]) * this->n_vy * this->n_vz +
                                    (j + CELL_VERTICES[l][1]) * this->n_vz +
                                    (k + CELL_VERTICES[l][2]);
                    vertex_signs[l]    = level_set[vertex_ids[l]] > 0;
                    vertices_coords[l] = {(CELL_VERTICES[l][0] + i) * cellWidth + origin.x(),
                                          (CELL_VERTICES[l][1] + j) * cellWidth + origin.y(),
                                          (CELL_VERTICES[l][2] + k) * cellWidth + origin.z()};
                }
                triangulation = get_marching_cubes_cell_triangulation(vertex_signs);
                for (int m = 0; m < 5; ++m) {
                    std::array<int, 3> triangle = triangulation[m];
                    if (triangle[0] >= 0) {
                        for (int l = 0; l < 3; ++l) {
                            vertexPair = CELL_EDGES[triangle[l]];
                            edgeIdx =
                                vertex_ids[vertexPair[0]] * 3 + CELL_EDGES_DIRECTION[triangle[l]];
                            alpha = level_set[vertex_ids[vertexPair[0]]] /
                                    (level_set[vertex_ids[vertexPair[0]]] -
                                     level_set[vertex_ids[vertexPair[1]]]);
                            intersection = (1.0 - alpha) * vertices_coords[vertexPair[0]] +
                                           alpha * vertices_coords[vertexPair[1]];
                            this->intersections[intersecId] = intersection;
                            this->edgeIntersection[edgeIdx] = intersecId;
                            this->triangles[triangleId][l]  = intersecId;
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
    this->intersections.resize(intersecId);
    this->triangles.resize(triangleId);
}

void learnSPH::surface::MarchingCubes::get_Isosurface_sparse(
    std::unordered_map<uint64_t, double> &level_map)
{
    std::unordered_set<uint64_t> cellIdx;
    std::array<uint, 8> vertexIds;
    std::array<bool, 8> vertexSigns;
    std::array<Eigen::Vector3d, 8> vertices_coords;
    std::array<std::array<int, 3>, 5> triangulation;
    std::array<int, 2> vertexPair;
    std::array<int64_t, 8> adjacentCells;
    std::array<int, 4> adjacentindices;
    std::vector<Eigen::Vector3d> gridVertices(this->n_vx * this->n_vy * this->n_vz);
    uint intersecId, triangleId, edgeIdx, va, vb;
    Eigen::Vector3d intersection;

    intersecId = 0;
    triangleId = 0;

    // create real coordinate mapping from vertex indices
    for (int i = 0; i < n_vx; i++) {
        for (int j = 0; j < n_vy; j++) {
            for (int k = 0; k < n_vz; k++) {
                uint vertexIndex = i * n_vy * n_vz + j * n_vz + k;
                gridVertices[vertexIndex] =
                    Eigen::Vector3d(i * cellWidth, j * cellWidth, k * cellWidth) + this->origin;
            }
        }
    }

    for (auto &idx : level_map) {
        if (idx.first > 4e9) {
            std::cout << "High Index found! Idx: " << idx.first << "\n";
            std::cout << level_map[idx.first] << "\n";
            exit(-1);
        }
        if (idx.second > 0) {
            for (int i = 0; i < 8; ++i) {
                adjacentCells[i] = learnSPH::utils::vertex8NeighborCells(idx.first, i, this->n_cx,
                                                                         this->n_cy, this->n_cz);
            }
            for (int i = 0; i < 6; i++) {
                int neighborIdx = learnSPH::utils::vertexSixNeighbors(
                    idx.first, i, gridVertices, this->n_vx, this->n_vy, this->n_vz);
                if (level_map.count(neighborIdx) > 0 && neighborIdx != -1) {
                    if (level_map[neighborIdx] > 0) {
                        // lookup for correct cell indices depending on neighbor index mapping
                        adjacentindices = learnSPH::utils::celladjByEdge(i);
                        for (int j = 0; j < 4; ++j) {
                            if (adjacentCells[adjacentindices[j]] != -1)
                                cellIdx.insert(adjacentCells[adjacentindices[j]]);
                        }
                    }
                }
            }
        }
    }

    for (auto &idx : cellIdx) {

        for (int i = 0; i < 8; i++) {
            vertexIds[i] = learnSPH::utils::cubeVertex2VertexIndex(idx, i, gridVertices, this->n_vx,
                                                                   this->n_vy, this->n_vz);
            if (vertexIds[i] != -1) {
                vertexSigns[i]     = level_map[vertexIds[i]] > 0;
                vertices_coords[i] = gridVertices[vertexIds[i]];
            }
        }
        triangulation = get_marching_cubes_cell_triangulation(vertexSigns);
        for (int i = 0; i < 5; i++) {
            std::array<int, 3> triangle = triangulation[i];
            if (triangle[0] >= 0) {
                for (int j = 0; j < 3; j++) {
                    // vertexPair = CELL_EDGES[triangle[j]];
                    va = CELL_EDGES[triangle[j]][0];
                    vb = CELL_EDGES[triangle[j]][1];
                    // edgeIdx    = vertexIds[vertexPair[0]] * 3 +
                    // CELL_EDGES_DIRECTION[triangle[j]]; double alpha =
                    //     level_map[vertexIds[vertexPair[0]]] /
                    //     (level_map[vertexIds[vertexPair[0]]] -
                    //     level_map[vertexIds[vertexPair[1]]]);

                    // Eigen::Vector3d intersection =
                    //     (1.0 - alpha) * gridVertices[vertexIds[vertexPair[0]]] +
                    //     alpha * vertices_coords[vertexPair[1]];
                    edgeIdx      = vertexIds[va] * 3 + CELL_EDGES_DIRECTION[triangle[j]];
                    double alpha = level_map[vertexIds[va]] /
                                   (level_map[vertexIds[va]] - level_map[vertexIds[vb]]);
                    intersection =
                        (1.0 - alpha) * vertices_coords[va] + alpha * vertices_coords[vb];
                    this->intersections[intersecId] = intersection;
                    this->triangles[triangleId][j]  = intersecId;
                    ++intersecId;
                }
                ++triangleId;
            } else {
                break;
            }
        }
    }
    this->intersections.resize(intersecId);
    this->triangles.resize(triangleId);
}

void learnSPH::surface::MarchingCubes::compute_normals(std::vector<Eigen::Vector3d> &positions,
                                                       std::vector<double> &densities,
                                                       Eigen::Vector3d bound_min)
{
    Eigen::Vector3d vertex;
    if (this->intersections.size() == 0) {
        std::cout << "No mesh vertices to compute normals for!"
                  << "\n";
        exit(-1);
    }

    this->intersectionNormals.resize(this->intersections.size());

    static Eigen::Vector3d ex = {1, 0, 0};
    static Eigen::Vector3d ey = {0, 1, 0};
    static Eigen::Vector3d ez = {0, 0, 1};

    for (int i = 0; i < this->intersections.size(); ++i) {
        this->intersectionNormals[i][0] =
            thetaFunction.singleSignedDistance(intersections[i] + epsilon * ex, positions,
                                               densities, bound_min) -
            thetaFunction.singleSignedDistance(intersections[i] - epsilon * ex, positions,
                                               densities, bound_min);
        this->intersectionNormals[i][1] =
            thetaFunction.singleSignedDistance(intersections[i] + epsilon * ey, positions,
                                               densities, bound_min) -
            thetaFunction.singleSignedDistance(intersections[i] - epsilon * ey, positions,
                                               densities, bound_min);
        this->intersectionNormals[i][2] =
            thetaFunction.singleSignedDistance(intersections[i] + epsilon * ez, positions,
                                               densities, bound_min) -
            thetaFunction.singleSignedDistance(intersections[i] - epsilon * ez, positions,
                                               densities, bound_min);

        this->intersectionNormals[i] = this->intersectionNormals[i] / 2.0 * epsilon;
        this->intersectionNormals[i].normalize();
    }
}

// void learnSPH::surface::MarchingCubes::compute_normals()
// {
//     Eigen::Vector3d vertex;
//     if (this->intersections.size() == 0) {
//         std::cout << "No mesh vertices to compute normals for!"
//                   << "\n";
//         exit(-1);
//     }

//     this->intersectionNormals.resize(this->intersections.size());

//     static Eigen::Vector3d ex = {1, 0, 0};
//     static Eigen::Vector3d ey = {0, 1, 0};
//     static Eigen::Vector3d ez = {0, 0, 1};

//     for (int i = 0; i < this->intersections.size(); ++i) {
//         this->intersectionNormals[i][0] =
//             thetaFunction.singleSignedDistance(intersections[i] + epsilon * ex) -
//             thetaFunction.singleSignedDistance(intersections[i] - epsilon * ex);
//         this->intersectionNormals[i][1] =
//             thetaFunction.singleSignedDistance(intersections[i] + epsilon * ey) -
//             thetaFunction.singleSignedDistance(intersections[i] - epsilon * ey);
//         this->intersectionNormals[i][2] =
//             thetaFunction.singleSignedDistance(intersections[i] + epsilon * ez) -
//             thetaFunction.singleSignedDistance(intersections[i] - epsilon * ez);

//         this->intersectionNormals[i] = this->intersectionNormals[i] / 2.0 * epsilon;
//         this->intersectionNormals[i].normalize();
//     }
// }
