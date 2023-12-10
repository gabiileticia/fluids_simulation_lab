#include <Eigen/Dense>
#include <string>
#ifndef TYPES
#define TYPES

namespace learnSPH
{
namespace types
{
struct boundary
{
    std::string filename;
    bool inner;
    Eigen::Vector3d min;
    Eigen::Vector3d max;
};
typedef double (*ImplicitSurface)(Eigen::Vector3d position, void* args);
} // namespace types
} // namespace learnSPH

#endif