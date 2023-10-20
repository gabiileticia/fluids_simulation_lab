#include "kernel.h"

double learnSPH::kernel::cubic_spline(const double q)
{
	assert(q >= 0.0);
	constexpr double alpha = 3.0 / (2.0 * PI);

	double result = 0;

	// branchless programming for the win
	result = (alpha * (2.0/3.0 - q * q + .5*q*q*q)) * (q < 1.0);
	result = (alpha * (1.0/6.0*(2 - q)*(2 - q)*(2 - q))) * ((1 <= q) && (q < 2.0));
	// result stays zero if q is greater equal than 2
	return result;
}

double learnSPH::kernel::cubic_grad_spline(const double q)
{
	// ...
	return 0.0;
}

