#ifndef MARCHING_CUBES
#define MARCHING_CUBES

#include <Eigen/Dense>
#include <array>
#include <cstdint>
#include <sys/types.h>
#include <unordered_map>
#include <vector>

#include "theta_functions.h"
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
    double c;
    bool implicitFlag; // false for < 0 inside and > 0 outside, true: for vice versa
    Eigen::Vector3d origin;
    
    std::vector<Eigen::Vector3d> intersections;
    std::vector<Eigen::Vector3d> intersectionNormals;
    std::unordered_map<uint, uint> edgeIntersection;
    std::vector<std::array<int, 3>> triangles;
    // for debuggin
    std::vector<Eigen::Vector3d> debug;

    // slearnSPH::types::ImplicitSurface implicitSurfaceFunction;
    learnSPH::theta_functions::FluidThetaFunction thetaFunction;
    void *funcArgs;

    MarchingCubes(double cellWidth, uint n_x, uint n_y, uint n_z, Eigen::Vector3d origin,
                  learnSPH::theta_functions::FluidThetaFunction &thetaFunction, double epsilon, bool implicitFlag, double c);

    void get_Isosurface_sparse(std::unordered_map<uint64_t, double> &level_map);
    void get_Isosurface(std::vector<double> &level_set);
    void get_Isosurface_for_real_this_time(std::vector<double> &level_set);
    void compute_normals(std::vector<Eigen::Vector3d> &positions, std::vector<double> &densities, Eigen::Vector3d bound_min);
    void compute_normals();
};
} // namespace surface
} // namespace learnSPH

#endif