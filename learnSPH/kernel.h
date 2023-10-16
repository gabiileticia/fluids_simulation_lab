#pragma once
#include <cassert>
#include <Eigen/Dense>

namespace learnSPH
{
	namespace kernel
	{
		constexpr double PI = 3.14159265358979323846;
		
		double cubic_spline(const double q);
		double kernel_function(const double xi, const double xj, const double h);
		double cubic_grad_spline(const double q);
	};
};