#include "kernel.h"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <math.h>

learnSPH::kernel::CubicSplineKernel::CubicSplineKernel(double h){
  h = h;
  h_inverse = 1 / h;
  gama = alpha / (h * h * h);
}

  double learnSPH::kernel::CubicSplineKernel::kernel_function(Eigen::Vector3d x){
  const double q = x.norm() * h_inverse;
    assert(q >= 0.0);
    double cubic = 0.0;

    if (q < 1.0) {
      cubic = (2.0 / 3.0 - q * q + 0.5 * q * q * q);
    } else if (q < 2.0) {
      cubic = (1.0 / 6.0) * (2 - q) * (2 - q) * (2 - q);
    }

    return gama * cubic;
}

Eigen::Vector3d learnSPH::kernel::CubicSplineKernel::kernel_gradient(Eigen::Vector3d x) {
    const double x_norm = x.norm();
    if (std::abs(x_norm - 0.0) < 0.000001) {
      return Eigen::Vector3d(0.0, 0.0, 0.0);
    }
    double q = x_norm * h_inverse;
    assert(q >= 0.0);

    double gama = alpha / (h * h * h * h * x_norm);
    double cubic = 0.0;

    if (q < 1.0) {
      cubic = (-2.0 * q + 1.5 * q * q);
    } else if (q < 2.0) {
      cubic = (-0.5) * (2 - q) * (2 - q);
    }

    return cubic * gama * x;
  }
