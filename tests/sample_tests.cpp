#include <cmath>
#include <cstdio>
#include <ostream>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <random>

#include "catch.hpp"

// uncomment for clangd to link the library
// #include "../extern/Eigen/Eigen/Dense"
#include "../learnSPH/kernel.h"
#include "../learnSPH/io.h"
#include "../learnSPH/sampling.h"
#include "../learnSPH/types.h"
#include "../learnSPH/time_integration.h"

#include <chrono>
#include <vector>


struct KernelFunction
{
    bool unity()
    {
        return true;
    }

    bool symmetry(learnSPH::kernel::CubicSplineKernel &kernel, Eigen::Vector3d xi, Eigen::Vector3d xj)
    {       
        return kernel.kernel_function(xi - xj) == kernel.kernel_function(xj - xi);
    }

    bool delta()
    {
        return true;
    }

    bool nonnegative(learnSPH::kernel::CubicSplineKernel &kernel, Eigen::Vector3d xi, Eigen::Vector3d xj)
    {
        return kernel.kernel_function(xi - xj) >= 0.0;
    }

    bool compactness(learnSPH::kernel::CubicSplineKernel &kernel, Eigen::Vector3d xi, Eigen::Vector3d xj, double beta, double tolerance)
    {
        if ((xi - xj).norm() > beta * kernel.h){
            return kernel.kernel_function(xi - xj) - 0 <= tolerance;
        }
        return true;
    }

    bool branch01_00(learnSPH::kernel::CubicSplineKernel &kernel,double tolerance)
    {
        // Branch 0-1
        // q = 0
        Eigen::Vector3d xi = Eigen::Vector3d(0.3, 0.3, 0.3);
        Eigen::Vector3d xj = Eigen::Vector3d(0.3, 0.3, 0.3);
        double q = (xi - xj).norm() / kernel.h;
        return std::abs(kernel.kernel_function(xi - xj) -0.0397887) < tolerance;
    }

    bool branch01_05(learnSPH::kernel::CubicSplineKernel &kernel,double tolerance)
    {
        // Branch 0-1
        // q = 0.5
        Eigen::Vector3d xi = Eigen::Vector3d(0.0, 0.0, 0.0);
        Eigen::Vector3d xj = Eigen::Vector3d(1.0, 0.0, 0.0);
        double q = (xi - xj).norm() / kernel.h;

        return std::abs(kernel.kernel_function(xi - xj) - 0.0285982) < tolerance;
    }

    bool branch12_10(learnSPH::kernel::CubicSplineKernel &kernel,double tolerance)
    {
        // Branch 1-2
        // q = 1.0
        Eigen::Vector3d xi = Eigen::Vector3d(-1.0, 0.0, 0.0);
        Eigen::Vector3d xj = Eigen::Vector3d(1.0, 0.0, 0.0);
        double q = (xi - xj).norm() / kernel.h;

        return std::abs(kernel.kernel_function(xi - xj) - 0.00994718) < tolerance;
    }
    
    bool branch12_15(learnSPH::kernel::CubicSplineKernel &kernel,double tolerance)
    {
        // Branch 1-2
        // q = 1.5
        Eigen::Vector3d xi = Eigen::Vector3d(-1.0, -1.0, 0.0);
        Eigen::Vector3d xj = Eigen::Vector3d(1.0, 1.0, 1.0);
        double q = (xi - xj).norm() / kernel.h;

        return std::abs(kernel.kernel_function(xi - xj) - 0.0012434) < tolerance;
    }

    bool branch2plus_20(learnSPH::kernel::CubicSplineKernel &kernel, double tolerance)
    {
        // Branch 2+
        // q = 2.0
        Eigen::Vector3d xi = Eigen::Vector3d(-1.0, 0.0, 0.0);
        Eigen::Vector3d xj = Eigen::Vector3d(3.0, 0.0, 0.0);
        double q = (xi - xj).norm() / kernel.h;

        return std::abs(kernel.kernel_function(xi - xj) - 0.0) < tolerance;
    }

    bool branch2plus_25(learnSPH::kernel::CubicSplineKernel &kernel, double tolerance)
    {
        // Branch 2+
        // q = 2.44949
        Eigen::Vector3d xi = Eigen::Vector3d(-1.0, -1.0, -1.0);
        Eigen::Vector3d xj = Eigen::Vector3d(3.0, 1.0, 1.0);
        double q = (xi - xj).norm() / kernel.h;

        return std::abs(kernel.kernel_function(xi - xj) - 0.0) < tolerance;
    }
};

struct GradientCubicSpline
{
    bool finite_differences_compare(learnSPH::kernel::CubicSplineKernel &kernel, Eigen::Vector3d xi, Eigen::Vector3d xj, double tolerance)
    {
        // Unit vectors
        Eigen::Vector3d ex = Eigen::Vector3d(1, 0, 0);
        Eigen::Vector3d ey = Eigen::Vector3d(0, 1, 0);
        Eigen::Vector3d ez = Eigen::Vector3d(0, 0, 1);

        //Finite differences individual components 
        double fin_diff_x = kernel.kernel_function(xi - xj + tolerance * ex) - kernel.kernel_function(xi - xj - tolerance * ex);
        double fin_diff_y = kernel.kernel_function(xi - xj + tolerance * ey) - kernel.kernel_function(xi - xj - tolerance * ey);
        double fin_diff_z = kernel.kernel_function(xi - xj + tolerance * ez) - kernel.kernel_function(xi - xj - tolerance * ez);

        //Finite differences
        Eigen::Vector3d fin_diff = Eigen::Vector3d(fin_diff_x, fin_diff_y, fin_diff_z) / (2.0 * tolerance) ;

        return (fin_diff - kernel.kernel_gradient(xi - xj)).norm() < tolerance;
    }

    bool gradient_zero_distance(learnSPH::kernel::CubicSplineKernel &kernel, Eigen::Vector3d x, double tolerance){
        return kernel.kernel_gradient(x - x) == Eigen::Vector3d(0, 0, 0);
    }
};

struct TimeIntegration {
    bool gravity_time_integration(learnSPH::timeIntegration::semiImplicitEuler &integrator, std::vector<learnSPH::types::Particle> &particles, double tolerance){
        return false;
    }
};

// Check out https://github.com/catchorg/Catch2 for more information about how to use Catch2
TEST_CASE( "Tests for our kernel function", "[kernel]" )
{
    KernelFunction kernel;
    GradientCubicSpline gradient;

    //Generate random number for h, and vectors xi, xj
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1, 1);

    double h = std::abs(dis(gen)) * 5;
    Eigen::Vector3d xi = Eigen::Vector3d(dis(gen), dis(gen), dis(gen));
    Eigen::Vector3d xj = Eigen::Vector3d(dis(gen), dis(gen), dis(gen));

    //Error variable 
    const double EPSILON = 0.000001;

    SECTION("Testing the branches of the kernel function") {

        learnSPH::kernel::CubicSplineKernel cubicSpline(2.0);

        REQUIRE(kernel.branch01_00(cubicSpline, EPSILON));
        REQUIRE(kernel.branch01_05(cubicSpline, EPSILON));
        REQUIRE(kernel.branch12_10(cubicSpline, EPSILON));
        REQUIRE(kernel.branch12_15(cubicSpline, EPSILON));
        REQUIRE(kernel.branch2plus_20(cubicSpline, EPSILON));
        REQUIRE(kernel.branch2plus_25(cubicSpline, EPSILON));
    }
    
    SECTION("Testing properties of the kernel function"){
        double beta = 2.0;
        double h = std::abs(dis(gen)) * 5;
        learnSPH::kernel::CubicSplineKernel cubicSpline2(h);

        for(int i = 0 ; i < 1000; i++)
        {
            Eigen::Vector3d xi = Eigen::Vector3d(dis(gen), dis(gen), dis(gen));
            Eigen::Vector3d xj = Eigen::Vector3d(dis(gen), dis(gen), dis(gen));
            
            REQUIRE(kernel.symmetry(cubicSpline2, xi,xj));
            REQUIRE(kernel.nonnegative(cubicSpline2, xi,xj));
            //doesn't work with class cubic spline
            REQUIRE(kernel.compactness(cubicSpline2, xi,xj,beta,EPSILON));
        }
        
        REQUIRE(kernel.unity());
        REQUIRE(kernel.delta());
    }

    SECTION("Testing gradient cubic spline") {
        for(int i = 0 ; i < 1000; i++)
        {
            double beta = 2.0;
            double h = std::abs(dis(gen)) * 5;
            Eigen::Vector3d xi = Eigen::Vector3d(dis(gen), dis(gen), dis(gen));
            Eigen::Vector3d xj = Eigen::Vector3d(dis(gen), dis(gen), dis(gen));

            learnSPH::kernel::CubicSplineKernel cubicSpline(h);

            REQUIRE(gradient.gradient_zero_distance(cubicSpline, xi, EPSILON));
            // Doesn't work with class cubicSpline
            REQUIRE(gradient.finite_differences_compare(cubicSpline, xi,xj, EPSILON));
        }
    }
}
TEST_CASE("Test for our time integration scheme. [semi-implicit Euler]"){
    double particle_radius = 0.25;
    // double particle_diameter = 2 * particle_radius;
    // double fluid_sampling_distance = particle_diameter; 
    // double boundary_sampling_distance = .8 * particle_diameter;
    // double smoothing_length = 1.2*particle_diameter;
    // double compact_support = 2.0 * smoothing_length;

    double dt = .5;
    double v_max = 981;

    double tolerance = 1e-7;

    TimeIntegration integrate_test;
    learnSPH::timeIntegration::semiImplicitEuler semImpEuler(dt, particle_radius, v_max);

    std::vector<learnSPH::types::Particle> particles(1000);

    SECTION("Testing time integration with gravity only."){
        REQUIRE(integrate_test.gravity_time_integration(semImpEuler,particles, tolerance))
    }

}