#include "kernel.h"
#include <math.h>
#include <iostream>

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

double learnSPH::kernel::kernel_function(const double xi, const double xj, const double h)
{
	const double q = fabs(xi - xj) / h;
	return (1 / (h*h*h)) * cubic_spline(q);
}

double learnSPH::kernel::cubic_grad_spline(const double q)
{
	// ...
	return 0.0;
}

