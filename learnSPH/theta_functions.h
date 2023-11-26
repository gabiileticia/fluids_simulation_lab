#include <Eigen/Dense>
#include "../learnSPH/kernel.h"

namespace learnSPH {
namespace theta_functions {
class FluidThetaFunction{
    public:
        double cell_width;
        double c;
        double support_radius;
        uint n_vy, n_vz;
        FluidThetaFunction(double c, double cell_width, double support_radius, uint n_vy, uint n_vz);
        // double ImplicitTorus(Eigen::Vector3d pos, void *args);
        void compute_fluid_reconstruction(std::vector<double> &level_set
                                        , std::vector<Eigen::Vector3d> positions, std::vector<double>  densities, Eigen::Vector3d bound_min
                                        , Eigen::Vector3d bound_max, learnSPH::kernel::CubicSplineKernel &kernel);
};
}
}