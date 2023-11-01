#include "kernel.h"
#include <cmath>
#include <math.h>
#include <iostream>
#include <Eigen/Dense>


double learnSPH::kernel::cubic_spline_branch(const double q)
{
	assert(q >= 0.0);
	constexpr double alpha = 3.0 / (2.0 * PI);
	double cubic = 0.0;

	if (q < 1.0) {
		cubic = alpha * (2.0/3.0 - q*q + 0.5 * q * q * q);
	}
	else if (q < 2.0) {
		cubic = alpha * (1.0/6.0) * (2 - q) * (2 - q) * (2 - q);
	}

	return alpha * cubic;
}

double learnSPH::kernel::cubic_spline(const double q){
	assert(q >= 0.0);
	constexpr double alpha = 3.0 / (2.0 * PI);

	double result = 0;

	result += (q < 1.0) * ((2.0 / 3.0) - q * q + 0.5 * q * q * q);
	result += ((q >= 1.0) && (q < 2.0)) * ((1.0 / 6.0) * (2 - q)*(2 - q)*(2 - q));

	return alpha * result;
}

double learnSPH::kernel::kernel_function(Eigen::Vector3d x, const double h)
{
	const double q = x.norm() / h;
	return (1 / (h*h*h)) * cubic_spline(q);
}

double learnSPH::kernel::kernel_function_with_cubic(Eigen::Vector3d x, const double h)
{
	const double q = x.norm() / h;
	assert(q >= 0.0);
	constexpr double alpha = 3.0 / (2.0 * PI);
	double gama = alpha / (h*h*h);
	double cubic = 0.0;

	if (q < 1.0) {
		cubic = (2.0/3.0 - q*q + 0.5*q*q*q);
	}
	else if (q < 2.0) {
		cubic = (1.0/6.0) * (2 - q) * (2 - q) * (2 - q);
	}

	return gama * cubic;
}

double learnSPH::kernel::kernel_function_with_cubic_no_branching(Eigen::Vector3d x, const double h)
{
	const double q = x.norm() / h;
	assert(q >= 0.0);
	constexpr double alpha = 3.0 / (2.0 * PI);
	double gama = alpha / (h*h*h);

	double result = 0.0;

	result += (q < 1.0) * ((2.0 / 3.0) - q * q + 0.5 * q * q * q);
	result += ((q >= 1.0) && (q < 2.0)) * ((1.0 / 6.0) * (2 - q)*(2 - q)*(2 - q));

	return gama * result;
}

double learnSPH::kernel::cubic_grad_spline(const double q)
{
	assert(q >= 0.0);
	constexpr double alpha = 3.0 / (2.0 * PI);

	double result = 0;

	result += (q < 1.0) * (- 2.0 * q + 1.5 * q * q);
	result += ((q >= 1.0) && (q < 2.0)) * (-0.5) * (2 - q) * (2 - q);

	return alpha * result;
}

Eigen::Vector3d learnSPH::kernel::kernel_gradient(Eigen::Vector3d x, const double h)
{
	const double x_norm = x.norm();
	if (std::abs(x_norm - 0.0) < 0.000001){
		return Eigen::Vector3d(0, 0, 0);
	}

	const double q = x_norm / h; 
	double del_w_del_q = (1 / (h*h*h)) * cubic_grad_spline(q);
	
	return x * del_w_del_q / (h * x_norm);
}

Eigen::Vector3d learnSPH::kernel::kernel_gradient_with_cubic(Eigen::Vector3d x, const double h)
{	
	const double x_norm = x.norm();
	if (std::abs(x_norm - 0.0) < 0.000001){
		return Eigen::Vector3d(0.0, 0.0, 0.0);
	}

	double q = x_norm / h;
	assert(q >= 0.0);

	constexpr double alpha = 3.0 / (2.0 * PI);
	double gama = alpha / (h*h*h*h*x_norm);
	double cubic = 0.0;

	if (q < 1.0) {
		cubic = (- 2.0 * q + 1.5 * q * q);
	}
	else if (q < 2.0) {
		cubic = (-0.5) * (2 - q) * (2 - q);
	}

	return cubic * gama * x;
}

Eigen::Vector3d learnSPH::kernel::kernel_gradient_with_cubic_nobranching(Eigen::Vector3d x, const double h)
{	
	const double x_norm = x.norm();
	if (std::abs(x_norm - 0.0) < 0.000001){
		return Eigen::Vector3d(0.0, 0.0, 0.0);
	}

	double q = x_norm / h;
	assert(q >= 0.0);

	constexpr double alpha = 3.0 / (2.0 * PI);
	double gama = alpha / (h * h * h * h *x_norm);
	double cubic = 0.0;

	cubic += (q < 1.0) * (- 2.0 * q + 1.5 * q * q);
	cubic += ((q >= 1.0) && (q < 2.0)) * (-0.5) * (2 - q) * (2 - q);

	return cubic * gama * x;
}