#include <Eigen/Dense>

namespace learnSPH {
    namespace theta_functions {
        double ImplicitTorus(Eigen::Vector3d pos, void *args);
        std::vector<double> fluid_reconstruction(std::vector<Eigen::Vector3d>  positions
                                    , std::vector<double>  densities, Eigen::Vector3d min, Eigen::Vector3d max);
    }
}