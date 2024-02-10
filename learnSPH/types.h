#include "../extern/Eigen/Eigen/Dense"
#include <string>
#ifndef TYPES
#define TYPES

namespace learnSPH
{
namespace types
{
struct object
{
    std::string filename;
    bool noCheck; // activate only if its a rectangular volume
    Eigen::Vector3d min;
    Eigen::Vector3d max;
    double lifetime = 0;
};
typedef double (*ImplicitSurface)(Eigen::Vector3d position, void *args);
struct emitter_data
{
    Eigen::Vector3d dir;
    Eigen::Vector3d origin;
    double r;
    double velocity;
    bool alternating;
    double emission_freq;
    int emit_counter;
};
struct fluid_sphere
{
    Eigen::Vector3d origin;
    double radius;
};
} // namespace types
} // namespace learnSPH

#endif