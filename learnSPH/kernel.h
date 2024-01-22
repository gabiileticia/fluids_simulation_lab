#pragma once
#include <cassert>
#include <Eigen/Dense>

#ifndef CUBIC_SPLINE
#define CUBIC_SPLINE

namespace learnSPH
{
	namespace kernel
	{
		constexpr double PI = 3.14159265358979323846;
		
		class CubicSplineKernel{
			public:
				static constexpr double alpha = 3.0 / (2.0 * learnSPH::kernel::PI);
				double gama;
				double rad_gamma;
				double h;
				double h_inverse;
				double c;
				CubicSplineKernel(double h, double beta);
				double kernel_function(Eigen::Vector3d x);
				Eigen::Vector3d kernel_gradient(Eigen::Vector3d x);
				double cohesion_kernel_function(double r);
				double adhesion_kernel_function(double r);
		};
	};
};

#endif