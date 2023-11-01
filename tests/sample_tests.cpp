#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <random>

#include "catch.hpp"

#include "../learnSPH/kernel.h"
#include "../learnSPH/io.h"
#include "../learnSPH/sampling.h"

#include <chrono>


struct KernelFunction
{
    bool unity()
    {
        return true;
    }

    bool symmetry(Eigen::Vector3d xi, Eigen::Vector3d xj, double h)
    {       
        return learnSPH::kernel::kernel_function(xi - xj,h) == learnSPH::kernel::kernel_function(xj - xi,h);
    }

    bool delta()
    {
        return true;
    }

    bool nonnegative(Eigen::Vector3d xi, Eigen::Vector3d xj, double h)
    {
        return learnSPH::kernel::kernel_function(xi - xj,h) >= 0.0;
    }

    bool compactness(Eigen::Vector3d xi, Eigen::Vector3d xj, double h, double beta, double tolerance)
    {
        if ((xi - xj).norm() > beta * h){
            return learnSPH::kernel::kernel_function(xi - xj,h) == 0.0;
        }
        return true;
    }

    bool branch01_00(double tolerance)
    {
        // Branch 0-1
        // q = 0
        Eigen::Vector3d xi = Eigen::Vector3d(0.3, 0.3, 0.3);
        Eigen::Vector3d xj = Eigen::Vector3d(0.3, 0.3, 0.3);
        double h = 2.0;
        double q = (xi - xj).norm() / h;
        return std::abs(learnSPH::kernel::kernel_function(xi - xj,h) -0.0397887) < tolerance;
    }

    bool branch01_05(double tolerance)
    {
        // Branch 0-1
        // q = 0.5
        Eigen::Vector3d xi = Eigen::Vector3d(0.0, 0.0, 0.0);
        Eigen::Vector3d xj = Eigen::Vector3d(1.0, 0.0, 0.0);
        double h = 2.0;
        double q = (xi - xj).norm() / h;

        return std::abs(learnSPH::kernel::kernel_function(xi - xj,h) - 0.0285982) < tolerance;
    }

    bool branch12_10(double tolerance)
    {
        // Branch 1-2
        // q = 1.0
        Eigen::Vector3d xi = Eigen::Vector3d(-1.0, 0.0, 0.0);
        Eigen::Vector3d xj = Eigen::Vector3d(1.0, 0.0, 0.0);
        double h = 2.0;
        double q = (xi - xj).norm() / h;

        return std::abs(learnSPH::kernel::kernel_function(xi - xj,h) - 0.00994718) < tolerance;
    }
    
    bool branch12_15(double tolerance)
    {
        // Branch 1-2
        // q = 1.5
        Eigen::Vector3d xi = Eigen::Vector3d(-1.0, -1.0, 0.0);
        Eigen::Vector3d xj = Eigen::Vector3d(1.0, 1.0, 1.0);
        double h = 2.0;
        double q = (xi - xj).norm() / h;

        return std::abs(learnSPH::kernel::kernel_function(xi - xj,h) - 0.0012434) < tolerance;
    }

    bool branch2plus_20(double tolerance)
    {
        // Branch 2+
        // q = 2.0
        Eigen::Vector3d xi = Eigen::Vector3d(-1.0, 0.0, 0.0);
        Eigen::Vector3d xj = Eigen::Vector3d(3.0, 0.0, 0.0);
        double h = 2.0;
        double q = (xi - xj).norm() / h;

        return std::abs(learnSPH::kernel::kernel_function(xi - xj,h) - 0.0) < tolerance;
    }

    bool branch2plus_25(double tolerance)
    {
        // Branch 2+
        // q = 2.44949
        Eigen::Vector3d xi = Eigen::Vector3d(-1.0, -1.0, -1.0);
        Eigen::Vector3d xj = Eigen::Vector3d(3.0, 1.0, 1.0);
        double h = 2.0;
        double q = (xi - xj).norm() / h;

        return std::abs(learnSPH::kernel::kernel_function(xi - xj,h) - 0.0) < tolerance;
    }
};

struct GradientCubicSpline
{
    bool finite_differences_compare(Eigen::Vector3d xi, Eigen::Vector3d xj, double h, double tolerance)
    {
        // Unit vectors
        Eigen::Vector3d ex = Eigen::Vector3d(1, 0, 0);
        Eigen::Vector3d ey = Eigen::Vector3d(0, 1, 0);
        Eigen::Vector3d ez = Eigen::Vector3d(0, 0, 1);

        //Finite differences individual components 
        double fin_diff_x = learnSPH::kernel::kernel_function(xi - xj + tolerance * ex,h) - learnSPH::kernel::kernel_function(xi - xj - tolerance * ex,h);
        double fin_diff_y = learnSPH::kernel::kernel_function(xi - xj + tolerance * ey,h) - learnSPH::kernel::kernel_function(xi - xj - tolerance * ey,h);
        double fin_diff_z = learnSPH::kernel::kernel_function(xi - xj + tolerance * ez,h) - learnSPH::kernel::kernel_function(xi - xj - tolerance * ez,h);

        //Finite differences
        Eigen::Vector3d fin_diff = Eigen::Vector3d(fin_diff_x, fin_diff_y, fin_diff_z) / (2.0 * tolerance) ;

        return (fin_diff - learnSPH::kernel::kernel_gradient(xi - xj,h)).norm() < tolerance;
    }

    bool gradient_zero_distance(Eigen::Vector3d x, double h, double tolerance){
        return learnSPH::kernel::kernel_gradient(x - x,h) == Eigen::Vector3d(0, 0, 0);
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
        REQUIRE(kernel.branch01_00(EPSILON));
        REQUIRE(kernel.branch01_05(EPSILON));
        REQUIRE(kernel.branch12_10(EPSILON));
        REQUIRE(kernel.branch12_15(EPSILON));
        REQUIRE(kernel.branch2plus_20(EPSILON));
        REQUIRE(kernel.branch2plus_25(EPSILON));
    }
    
    SECTION("Testing properties of the kernel function"){
        double beta = 2.0;
        for(int i = 0 ; i < 1000; i++)
        {
            double h = std::abs(dis(gen)) * 5;
            Eigen::Vector3d xi = Eigen::Vector3d(dis(gen), dis(gen), dis(gen));
            Eigen::Vector3d xj = Eigen::Vector3d(dis(gen), dis(gen), dis(gen));

            REQUIRE(kernel.compactness(xi,xj,h,beta,EPSILON));
            REQUIRE(kernel.symmetry(xi,xj,h));
            REQUIRE(kernel.nonnegative(xi,xj,h));
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

            REQUIRE(gradient.finite_differences_compare(xi,xj,h,EPSILON));
            REQUIRE(gradient.gradient_zero_distance(xi,h,EPSILON));
        }
    }
}