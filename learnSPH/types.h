#include <Eigen/Dense>
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
};
typedef double (*ImplicitSurface)(Eigen::Vector3d position, void* args);
} // namespace types
} // namespace learnSPH

#endif