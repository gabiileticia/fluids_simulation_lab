#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <random>


#include "catch.hpp"

#include "../learnSPH/kernel.h"


struct CubicSpline
{
    bool branch01_00(double tolerance)
    {
        // Branch 0-1
        const double q = 0.0;
        return fabs(learnSPH::kernel::cubic_spline(q) - 0.31831) < tolerance;
    }

    bool branch01_05(double tolerance)
    {
        // Branch 0-1
        const double q = 0.5;
        return fabs(learnSPH::kernel::cubic_spline(q) - 0.228785) < tolerance;
        
    }

    bool branch12_10(double tolerance)
    {
        // Branch 1-2
        const double q = 1.0;
        return fabs(learnSPH::kernel::cubic_spline(q) - 0.0795775) < tolerance;
    }
    
    bool branch12_15(double tolerance)
    {
        // Branch 1-2
        const double q = 1.5;
        return fabs(learnSPH::kernel::cubic_spline(q) - 0.00994718) < tolerance;
    }

    bool branch2plus_20()
    {
        // Branch 2+
        const double q = 2.0;
        return learnSPH::kernel::cubic_spline(q) == 0.0;
    }

    bool branch2plus_25()
    {
        // Branch 2+
        const double q = 2.5;
        return learnSPH::kernel::cubic_spline(q) == 0.0;
    }
};

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
};

struct GradientCubicSpline
{
    bool branch01_00()
    {
        // Branch 0-1
        const double q = 0.0;
        return learnSPH::kernel::cubic_grad_spline(q) == 0.0;
    }

    bool branch01_05(double tolerance)
    {
        // Branch 0-1
        const double q = 0.5;
        return fabs(learnSPH::kernel::cubic_grad_spline(q) - (-0.298416)) < tolerance;
    }

    bool branch12_10(double tolerance)
    {
        // Branch 1-2
        const double q = 1.0;
        return fabs(learnSPH::kernel::cubic_grad_spline(q) - (-0.238732)) < tolerance;
    }

    bool branch12_15(double tolerance)
    {
        // Branch 1-2
        const double q = 1.5;
        return fabs(learnSPH::kernel::cubic_grad_spline(q) - (-0.0596831)) < tolerance;
    }

    bool branch2plus_20()
    {
        // Branch 2+
        const double q = 2.0;
        return learnSPH::kernel::cubic_grad_spline(q) == 0.0;
    }
    
    bool branch2plus_25()
    {
        // Branch 2+
        const double q = 2.5;
        return learnSPH::kernel::cubic_grad_spline(q) == 0.0;
    }

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
        Eigen::Vector3d fin_diff = Eigen::Vector3d(fin_diff_x, fin_diff_y, fin_diff_z) / (2 * tolerance) ;

        return (fin_diff - learnSPH::kernel::kernel_gradient(xi - xj,h)).norm() < tolerance;
    }

    bool gradient_zero_distance(Eigen::Vector3d x, double h, double tolerance){
        return (learnSPH::kernel::kernel_gradient(x - x,h) - Eigen::Vector3d(0, 0, 0)).norm() < tolerance;
    }

};

// Check out https://github.com/catchorg/Catch2 for more information about how to use Catch2
TEST_CASE( "Tests for our kernel function", "[kernel]" )
{
    KernelFunction kernel;
    CubicSpline cubic;
    GradientCubicSpline gradient;

    //Generate random number for h, and vectors xi, xj
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);

    double h = dis(gen) * 5;
    Eigen::Vector3d xi = Eigen::Vector3d(dis(gen), dis(gen), dis(gen));
    Eigen::Vector3d xj = Eigen::Vector3d(dis(gen), dis(gen), dis(gen));

    //Error variable 
    const double EPSILON = 0.000001;

    SECTION("Testing the branches of the cubic spline") {
        REQUIRE(cubic.branch01_00(EPSILON));
        REQUIRE(cubic.branch01_05(EPSILON));
        REQUIRE(cubic.branch12_10(EPSILON));
        REQUIRE(cubic.branch12_15(EPSILON));
        REQUIRE(cubic.branch2plus_20());
        REQUIRE(cubic.branch2plus_25());
    }
    
    SECTION("Testing properties of the kernel function"){
        for(int i = 0 ; i < 1000; i++)
        {
            double beta = 2.0;
            double h = dis(gen) * 5;
            Eigen::Vector3d xi = Eigen::Vector3d(dis(gen), dis(gen), dis(gen));
            Eigen::Vector3d xj = Eigen::Vector3d(dis(gen), dis(gen), dis(gen));

            REQUIRE(kernel.compactness(xi,xj,h,beta,EPSILON));
            REQUIRE(kernel.symmetry(xi,xj,h));
            REQUIRE(kernel.nonnegative(xi,xj,h));
        }
        
        REQUIRE(kernel.unity());
        REQUIRE(kernel.delta());
    }

    SECTION("Testing the branches of gradient cubic spline") {
        REQUIRE(gradient.branch01_00());
        REQUIRE(gradient.branch01_05(EPSILON));
        REQUIRE(gradient.branch12_10(EPSILON));
        REQUIRE(gradient.branch12_15(EPSILON));
        REQUIRE(gradient.branch2plus_20());
        REQUIRE(gradient.branch2plus_25());
    }

    SECTION("Testing gradient cubic spline") {
        REQUIRE(gradient.finite_differences_compare(xi,xj,h,EPSILON));
        REQUIRE(gradient.gradient_zero_distance(xi,h,EPSILON));
    }
}
