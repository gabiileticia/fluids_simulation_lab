#include "kernel.h"
#include <cmath>
#include <math.h>
#include <iostream>
#include <Eigen/Dense>

double learnSPH::kernel::cubic_spline(const double q)
{
	assert(q >= 0.0);
	std::cout << "\n\nQ value: " << q;
	constexpr double alpha = 3.0 / (2.0 * PI);

	double result = 0;
	// double branch1 = q < 1.0;
	// double branch2 = (1 <= q)&&(q < 2.0);

	// double branch1_result = (alpha * ((2.0/3.0) - q * q + .5*std::pow(q,3)));
	// double branch2_result = (alpha * ((1.0/6.0)*std::pow(2-q, 3)));

	// // branchless programming for the win
	// result = branch1 * branch1_result;
	// std::cout << "\nBranch 1\n\tbool: " << branch1;
	// std::cout << "\n\tresult: " << branch1_result;
	// std::cout << "\n\tresults multiplied with bool: " << result;

	// result = branch2 * branch2_result;
	// std::cout << "\nBranch 2\n\tbool: " << branch2;
	// std::cout << "\n\tresult: " << branch2_result;
	// std::cout << "\n\tresults multiplied with bool: " << result;
	// result stays zero if q is greater equal than 2

	if (q < 1.0) {
		result = alpha * ((2.0 / 3.0) - std::pow(q, 2) * .5 * std::pow(q, 3));
	} else if (q < 2.0) {
		result = alpha * ((1.0 / 6.0) * std::pow(2 - q, 3));
	} else {
		result = 0;
	}
	return result;
}

double learnSPH::kernel::kernel_function(Eigen::Vector3d x, const double h)
{
	const double q = x.norm() / h;
	return (1 / (std::pow(h, 3))) * cubic_spline(q);
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
	const double x_norm = x.norm();
	if (x_norm == 0.0){
		return Eigen::Vector3d(0, 0, 0);
	}

	const double q = x_norm / h; 
	double del_w_del_q = (1 / (h*h*h)) * cubic_grad_spline(q);
	
	return x * del_w_del_q / (h * x_norm);
}

