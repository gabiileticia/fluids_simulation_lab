#include <array>
#include <chrono>
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

learnSPH::surface::MarchingCubes::MarchingCubes(double cellWidth, uint n_x, uint n_y, uint n_z,
                                                Eigen::Vector3d origin, double epsilon)
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

    this->triangles.resize(n_vx * n_vy * n_vz * 5);
    this->intersections.resize(this->triangles.size() * 3);
}

void learnSPH::surface::MarchingCubes::get_Isosurface(std::vector<double> &level_set)
{
    uint vertex_id_a, vertex_id_b, edge_id, intersecId, triangleId;
    std::array<uint, 8> vertex_ids;
    Eigen::Vector3d intersection;
    std::array<bool, 8> vertex_signs;
    std::array<std::array<int, 3>, 5> triangulation;
    std::array<int, 2> vertex_pair;
    double alpha;
    Eigen::Vector3d vertex_coord_a;
    Eigen::Vector3d vertex_coord_b;

    intersecId = 0;
    triangleId = 0;

    std::vector<Eigen::Vector3d> directions(3);
    directions[0] = Eigen::Vector3d(1, 0, 0);
    directions[1] = Eigen::Vector3d(0, 1, 0);
    directions[2] = Eigen::Vector3d(0, 0, 1);

    for (int i = 0; i < n_vx; i++) {
        for (int j = 0; j < n_vy; j++) {
            for (int k = 0; k < n_vz; k++) {

                vertex_id_a = i * this->n_vy * this->n_vz + j * this->n_vz + k;
                vertex_coord_a =
                    Eigen::Vector3d(i * cellWidth, j * cellWidth, k * cellWidth) + origin;

                for (int l = 0; l < 3; l++) {

                    vertex_id_b = (i + directions[l].x()) * this->n_vy * this->n_vz +
                                  (j + directions[l].y()) * this->n_vz + (k + directions[l].z());
                    vertex_coord_b = Eigen::Vector3d((i + directions[l].x()) * cellWidth,
                                                     (j + directions[l].y()) * cellWidth,
                                                     (k + directions[l].z()) * cellWidth) +
                                     origin;

                    if ((vertex_id_b > level_set.size() - 1) || vertex_id_b < 0) {
                        continue;
                    }

                    if ((level_set[vertex_id_a] > epsilon && level_set[vertex_id_b] < -epsilon) ||
                        (level_set[vertex_id_a] < -epsilon && level_set[vertex_id_b] > epsilon)) {

                        edge_id = vertex_id_a * 3 + l;
                        alpha   = level_set[vertex_id_a] /
                                (level_set[vertex_id_a] - level_set[vertex_id_b]);
                        intersection = (1.0 - alpha) * vertex_coord_a + alpha * vertex_coord_b;

                        this->intersections[intersecId] = intersection;
                        this->edgeIntersection[edge_id] = intersecId;
                        ++intersecId;
                    }
                }
            }
        }
    }

    for (int i = 0; i < this->n_cx; i++) {
        for (int j = 0; j < this->n_cy; j++) {
            for (int k = 0; k < this->n_cz; k++) {

                for (int l = 0; l < 8; l++) {
                    vertex_ids[l] = (i + CELL_VERTICES[l][0]) * this->n_vy * this->n_vz +
                                    (j + CELL_VERTICES[l][1]) * this->n_vz +
                                    (k + CELL_VERTICES[l][2]);
                    vertex_signs[l] = level_set[vertex_ids[l]] < 0;
                }

                triangulation = get_marching_cubes_cell_triangulation(vertex_signs);

                for (int m = 0; m < 5; m++) {
                    std::array<int, 3> triangle = triangulation[m];
                    if (triangle[0] >= 0) {
                        for (int l = 0; l < 3; ++l) {
                            vertex_pair = CELL_EDGES[triangle[l]];
                            edge_id =
                                vertex_ids[vertex_pair[0]] * 3 + CELL_EDGES_DIRECTION[triangle[l]];
                            this->triangles[triangleId][l] = edgeIntersection[edge_id];
                        }
                        triangleId++;
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
    uint x, y, z;
    std::unordered_set<uint64_t> cellIdx;
    std::array<uint, 8> vertexIds;
    std::array<bool, 8> vertexSigns;
    std::array<Eigen::Vector3d, 8> vertices_coords;
    std::array<std::array<int, 3>, 5> triangulation;
    std::array<int, 2> vertexPair;
    std::array<int64_t, 8> adjacentCells;
    std::array<int, 4> adjacentindices;
    uint intersecId, triangleId, edgeIdx;
    Eigen::Vector3d intersection, vertexCoord, va, vb;
    double alpha;

    intersecId = 0;
    triangleId = 0;

    debug.clear();

    for (auto &idx : level_map) {
        if (idx.first > 4e9) {
            std::cout << "To high Index found! Idx: " << idx.first << "\n";
            std::cout << level_map[idx.first] << "\n";
            exit(-1);
        }

        if (idx.second < 0) {
            for (int i = 0; i < 8; ++i) {
                adjacentCells[i] = learnSPH::utils::vertex8NeighborCells(idx.first, i, this->n_vx,
                                                                         this->n_vy, this->n_vz);
            }
            for (int i = 0; i < 6; i++) {
                int neighborIdx = learnSPH::utils::vertexSixNeighbors(idx.first, i, this->n_vx,
                                                                      this->n_vy, this->n_vz);
                adjacentindices = learnSPH::utils::celladjByEdge(i);
                if (neighborIdx != -1) {
                    if (level_map.count(neighborIdx) > 0) {
                        if (level_map[neighborIdx] > 0) {
                            // lookup for correct cell indices depending on neighbor index mapping
                            // computing intersection positions
                            if (i < 3) {
                                edgeIdx = idx.first * 3 + i % 3;
                                alpha   = level_map[idx.first] /
                                        (level_map[idx.first] - level_map[neighborIdx]);
                                va = learnSPH::utils::index2coord(idx.first, cellWidth, this->n_vx,
                                                                  this->n_vy, this->n_vz,
                                                                  this->origin);
                                vb = learnSPH::utils::index2coord(neighborIdx, cellWidth,
                                                                  this->n_vx, this->n_vy,
                                                                  this->n_vz, this->origin);
                            } else {
                                edgeIdx = neighborIdx * 3 + i % 3;
                                alpha   = level_map[neighborIdx] /
                                        (level_map[neighborIdx] - level_map[idx.first]);
                                va = learnSPH::utils::index2coord(neighborIdx, cellWidth,
                                                                  this->n_vx, this->n_vy,
                                                                  this->n_vz, this->origin);
                                vb = learnSPH::utils::index2coord(idx.first, cellWidth, this->n_vx,
                                                                  this->n_vy, this->n_vz,
                                                                  this->origin);
                            }
                            intersection                    = (1.0 - alpha) * va + alpha * vb;
                            this->intersections[intersecId] = intersection;
                            this->edgeIntersection[edgeIdx] = intersecId;
                            intersecId++;
                            // saving relevant cells for triangulation
                            for (int j = 0; j < 4; ++j) {
                                if (adjacentCells[adjacentindices[j]] != -1) {
                                    cellIdx.insert(adjacentCells[adjacentindices[j]]);
                                }
                            }
                        }
                    }
                } else {
                    for (int j = 0; j < 4; j++) {
                        adjacentCells[adjacentindices[j]] = -1;
                    }
                }
            }
        }
    }

    // marching through the cubes for triangulation
    for (auto &idx : cellIdx) {

        for (int i = 0; i < 8; i++) {
            vertexIds[i] =
                learnSPH::utils::cubeVertex2VertexIndex(idx, i, this->n_vx, this->n_vy, this->n_vz);
            if (vertexIds[i] != -1) {
                vertexSigns[i] = level_map[vertexIds[i]] < 0;
            }
        }

        triangulation = get_marching_cubes_cell_triangulation(vertexSigns);
        for (int i = 0; i < 5; i++) {
            std::array<int, 3> triangle = triangulation[i];
            if (triangle[0] >= 0) {
                for (int j = 0; j < 3; j++) {
                    vertexPair = CELL_EDGES[triangle[j]];
                    edgeIdx    = vertexIds[vertexPair[0]] * 3 + CELL_EDGES_DIRECTION[triangle[j]];
                    this->triangles[triangleId][j] = edgeIntersection[edgeIdx];
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

void learnSPH::surface::MarchingCubes::compute_normals()
{

    if (this->intersections.size() == 0) {
        std::cout << "No mesh vertices to compute normals for!"
                  << "\n";
    }

    Eigen::Vector3d normal;

    this->intersectionNormals.resize(this->intersections.size(), Eigen::Vector3d(0, 0, 0));
    std::vector<int> sum_triangles(intersectionNormals.size(), 0);

    for (int i = 0; i < triangles.size(); i++) {
        for (int j = 0; j < 3; j++) {

            Eigen::Vector3d u =
                intersections[triangles[i][(j + 1) % 3]] - intersections[triangles[i][j]];
            Eigen::Vector3d v =
                intersections[triangles[i][(j + 2) % 3]] - intersections[triangles[i][j]];

            normal[0] = u.y() * v.z() - u.z() * v.y();
            normal[1] = u.z() * v.x() - u.x() * v.z();
            normal[2] = u.x() * v.y() - u.y() * v.x();

            this->intersectionNormals[triangles[i][j]] += normal.normalized();
            sum_triangles[triangles[i][j]] += 1;
        }
    }
    for (int i = 0; i < intersectionNormals.size(); i++) {
        this->intersectionNormals[i] =
            (this->intersectionNormals[i] / sum_triangles[i]).normalized();
    }
}
