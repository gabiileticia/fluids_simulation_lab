#include "kernel.h"
#include <math.h>
#include <iostream>
#include <Eigen/Dense>

double learnSPH::kernel::cubic_spline(const double q)
{
	assert(q >= 0.0);
	constexpr double alpha = 3.0 / (2.0 * PI);

	if (q < 1.0) {
		return alpha * (2.0/3.0 - q*q + 0.5*q*q*q);
	}
	else if (q < 2.0) {
		return alpha * (1.0/6.0) * (2 - q) * (2 - q) * (2 - q);
	}
	else {
		return 0.0;
	}
}

double learnSPH::kernel::kernel_function(Eigen::Vector3d x, const double h)
{
	const double q = x.norm() / h;
	return (1 / (h*h*h)) * cubic_spline(q);
}

double learnSPH::kernel::cubic_grad_spline(const double q)
{
	assert(q >= 0.0);
	constexpr double alpha = 3.0 / (2.0 * PI);

	if (q < 1.0) {
		return alpha * (- 2.0*q + 1.5*q*q);
	}
	else if (q < 2.0) {
		return alpha * (-0.5) * (2 - q) * (2 - q);
	}
	else {
		return 0.0;
	}
}

Eigen::Vector3d learnSPH::kernel::kernel_gradient(Eigen::Vector3d x, const double h)
{
	const double q = x.norm() / h; // Maybe pass the norm already computed to the function would be better
	double del_w_del_q = (1 / (h*h*h)) * cubic_grad_spline(q);
	if (x.norm() == 0.0){
		return Eigen::Vector3d(0, 0, 0);
	}
	return del_w_del_q * x / (h * x.norm());
}

