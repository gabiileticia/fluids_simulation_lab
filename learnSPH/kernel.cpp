#include "kernel.h"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <math.h>

learnSPH::kernel::CubicSplineKernel::CubicSplineKernel(double h, double beta){
  this->h = h;
  this->h_inverse = 1 / h;
  this->c = beta;
  this->gama = alpha / (h * h * h);
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

double learnSPH::kernel::CubicSplineKernel::cohesion_kernel_function(double r){

    assert(r >= 0.0);
    double aux = 32.0/(learnSPH::kernel::PI * std::pow(this->c, 9));
    double kernel_branch = 0.0;
  
    if (r < (this->c / 2.0)) {
      kernel_branch = 2.0 * std::pow((this->c - r), 3) * std::pow(r, 3) - std::pow(c, 6) / 64;
    } else if (r < this->c) {
      kernel_branch = std::pow((this->c - r), 3) * std::pow(r,3);
    }
    return aux * kernel_branch;
}


double learnSPH::kernel::CubicSplineKernel::adhesion_kernel_function(double r){

    double aux = 0.007/std::pow(this->c, 3.25);
    double kernel_branch = 0.0;
  
    if (r >= (this->c / 2.0) && r<= this->c) {
      kernel_branch = std::pow(- 4.0 * r * r / c + 6.0 * r - 2.0 * c, 1.0/4.0);
    }

    return aux * kernel_branch;
}
