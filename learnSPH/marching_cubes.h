#ifndef MARCHING_CUBES
#define MARCHING_CUBES

#include <Eigen/Dense>
#include <array>
#include <sys/types.h>
#include <unordered_map>
#include <vector>

#include "types.h"

namespace learnSPH
{
namespace surface
{
class MarchingCubes
{
  public:
    uint n_x, n_z, n_y, n_vx, n_vy, n_vz, n_cx, n_cy, n_cz, n_ex, n_ey, n_ez;
    double cellWidth;
    double epsilon;
    bool implicitFlag; // false for < 0 inside and > 0 outside, true: for vice versa
    Eigen::Vector3d origin;

    // std::vector<double> levelSet; // Stores the SD value for position with value btween -1 (inside
    //                               // surface) and 1 (outside surface)
    std::vector<Eigen::Vector3d> intersections;
    std::vector<Eigen::Vector3d> intersectionNormals;
    std::unordered_map<uint, uint> edgeIntersection;
    std::vector<std::array<int, 3UL>> triangles;
    std::vector<std::array<int, 3>> triangles2;

    learnSPH::types::ImplicitSurface implicitSurfaceFunction;
    void *funcArgs;

    // MarchingCubes(double cellWidth, uint n_x, uint n_y, uint n_z, Eigen::Vector3d origin,
    //               learnSPH::types::ImplicitSurface implicitSurfaceFunction, void *funcArgs,
    //               double epsilon, bool implicitFlag);

    MarchingCubes(
    double cellWidth, uint n_x, uint n_y, uint n_z, Eigen::Vector3d origin, double epsilon, bool implicitFlag);

    void get_isosurface(std::vector<double> &level_set);
    void compute_normals();
    // void compute_normals_alternative();
    void compute_normals_gl(std::vector<double> &level_set);

};
} // namespace surface
} // namespace learnSPH

#endif