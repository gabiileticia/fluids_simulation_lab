#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <random>


#include "catch.hpp"

#include "../learnSPH/kernel.h"
#include "../learnSPH/neighborhood_search.h"
#include "../learnSPH/io.h"
#include "../learnSPH/sampling.h"
#include "../extern/CompactNSearch/include/CompactNSearch/CompactNSearch"


#include <chrono>


struct NeighborhoodSearch
{
    bool branch01_00(double h, std::vector<Eigen::Vector3d> particles, double beta)
    {
        auto start_brute = std::chrono::high_resolution_clock::now();

        learnSPH::neighborhood_search::brute_force_search(h, particles, beta);

        auto stop_brute = std::chrono::high_resolution_clock::now();
        auto duration_brute = std::chrono::duration_cast<std::chrono::microseconds>(stop_brute - start_brute);
        std::cout << "Brute force time: ";
        std::cout << duration_brute.count();
        std::cout << " microseconds." << std::endl;


        CompactNSearch::NeighborhoodSearch nsearch(beta);
        unsigned int point_set_id = nsearch.add_point_set(particles.front().data(), particles.size());

        auto start_compact = std::chrono::high_resolution_clock::now();
        nsearch.find_neighbors();

        auto stop_compact = std::chrono::high_resolution_clock::now();
        auto duration_compact = std::chrono::duration_cast<std::chrono::microseconds>(stop_compact - start_compact);
        std::cout << "CompactNSearch time: ";
        std::cout << duration_compact.count();
        std::cout << " microseconds." << std::endl;
        return true;
    }
};


// TEST_CASE( "Tests for neighborhood search", "[neighborhood search]" )
// {
//     NeighborhoodSearch neighbor;

//     // Load a obj surface mesh
// 	const std::vector<learnSPH::TriMesh> meshes = learnSPH::read_tri_meshes_from_obj("./res/box.obj");
// 	const learnSPH::TriMesh& box = meshes[0];

//     double h = 2.0;
//     double beta = 0.1;

//     SECTION("Testing the branches of the cubic spline") {
//         // Sample the mesh with particles
//         double sampling_distance = 0.5;
//         std::vector<Eigen::Vector3d> particles;

//         for(int i = 0 ; i < 5; i++){
//             learnSPH::sampling::triangle_mesh(particles, box.vertices, box.triangles, sampling_distance);

//             std::cout << "Number of particles: ";
// 	        std::cout << particles.size() << std::endl;

//             REQUIRE(neighbor.branch01_00(h, particles, beta));
//             sampling_distance = sampling_distance/2;
//         }
//     }
// }



struct CubicSpline
{
    bool improvement(int num_ite, std::mt19937 gen, double tolerance)
    {
        std::uniform_real_distribution<> dis_q(0, 5);
        auto start_nobranch = std::chrono::high_resolution_clock::now();
    
        for(int i = 0 ; i < num_ite; i++)
        {
            double q = dis_q(gen);
            learnSPH::kernel::cubic_spline(q);
        }

        auto stop_nobranch = std::chrono::high_resolution_clock::now();
        auto duration_nobranch = std::chrono::duration_cast<std::chrono::microseconds>(stop_nobranch - start_nobranch);
        std::cout << "Cubic Spline no branching time: ";
        std::cout << duration_nobranch.count();
        std::cout << " microseconds." << std::endl;


        auto start_branch = std::chrono::high_resolution_clock::now();

        for(int i = 0 ; i < num_ite; i++)
        {
            double q = dis_q(gen);
            learnSPH::kernel::cubic_spline_branch(q);
        }

        auto stop_branch = std::chrono::high_resolution_clock::now();
        auto duration_branch = std::chrono::duration_cast<std::chrono::microseconds>(stop_branch - start_branch);
        std::cout << "Cubic Spline with branching time: ";
        std::cout << duration_branch.count();
        std::cout << " microseconds." << std::endl;

        // return duration_nobranch.count() < duration_branch.count();
        return true;
    }
};


struct KernelFunction
{
    bool compare_results(Eigen::Vector3d xi, Eigen::Vector3d xj, double h, double tolerance){
        double k1 = learnSPH::kernel::kernel_function(xi - xj,h);
        double k2 = learnSPH::kernel::kernel_function_with_cubic(xi - xj,h);
        double k3 = learnSPH::kernel::kernel_function_with_cubic_no_branching(xi - xj,h);

        return std::abs(k1 - k2) < tolerance && std::abs(k1 - k3) < tolerance && std::abs(k2 - k3) < tolerance;
    }

    bool improvement(int num_ite, std::mt19937 gen, double tolerance)
    {
        std::uniform_real_distribution<> dis_k(-1, 1);

        auto start_simplekernel = std::chrono::high_resolution_clock::now();
    
        for(int i = 0 ; i < num_ite; i++)
        {
            double h = std::abs(dis_k(gen)) * 5;
            Eigen::Vector3d xi = Eigen::Vector3d(dis_k(gen), dis_k(gen), dis_k(gen));
            Eigen::Vector3d xj = Eigen::Vector3d(dis_k(gen), dis_k(gen), dis_k(gen));
            learnSPH::kernel::kernel_function(xi - xj,h);
            
        }

        auto stop_simplekernel = std::chrono::high_resolution_clock::now();
        auto duration_simplekernel = std::chrono::duration_cast<std::chrono::microseconds>(stop_simplekernel - start_simplekernel);
        std::cout << "Simple Kernel time: ";
        std::cout << duration_simplekernel.count();
        std::cout << " microseconds." << std::endl;


        auto start_nofunction = std::chrono::high_resolution_clock::now();

        for(int i = 0 ; i < num_ite; i++)
        {
            double h = std::abs(dis_k(gen)) * 5;
            Eigen::Vector3d xi = Eigen::Vector3d(dis_k(gen), dis_k(gen), dis_k(gen));
            Eigen::Vector3d xj = Eigen::Vector3d(dis_k(gen), dis_k(gen), dis_k(gen));
            learnSPH::kernel::kernel_function_with_cubic(xi - xj,h);
            
        }

        auto stop_nofunction = std::chrono::high_resolution_clock::now();
        auto duration_nofunction = std::chrono::duration_cast<std::chrono::microseconds>(stop_nofunction - start_nofunction);
        std::cout << "No function call kernel time: ";
        std::cout << duration_nofunction.count();
        std::cout << " microseconds." << std::endl;


        auto start_nobranch = std::chrono::high_resolution_clock::now();

        for(int i = 0 ; i < num_ite; i++)
        {
            double h = std::abs(dis_k(gen)) * 5;
            Eigen::Vector3d xi = Eigen::Vector3d(dis_k(gen), dis_k(gen), dis_k(gen));
            Eigen::Vector3d xj = Eigen::Vector3d(dis_k(gen), dis_k(gen), dis_k(gen));
            learnSPH::kernel::kernel_function_with_cubic_no_branching(xi - xj,h);
            
        }

        auto stop_nobranch = std::chrono::high_resolution_clock::now();
        auto duration_nobranch = std::chrono::duration_cast<std::chrono::microseconds>(stop_nobranch - start_nobranch);
        std::cout << "No function call and no branching kernel time: ";
        std::cout << duration_nobranch.count();
        std::cout << " microseconds." << std::endl;

        return true;
    }
};

struct GradientCubicSpline
{
    bool improvement(int num_ite, std::mt19937 gen, double tolerance)
    {

        std::uniform_real_distribution<> dis_k(-1, 1);

        auto start_simplekernel = std::chrono::high_resolution_clock::now();
    
        for(int i = 0 ; i < num_ite; i++)
        {
            double h = std::abs(dis_k(gen)) * 5;
            Eigen::Vector3d xi = Eigen::Vector3d(dis_k(gen), dis_k(gen), dis_k(gen));
            Eigen::Vector3d xj = Eigen::Vector3d(dis_k(gen), dis_k(gen), dis_k(gen));
            learnSPH::kernel::kernel_gradient(xi - xj,h);
            
        }

        auto stop_simplekernel = std::chrono::high_resolution_clock::now();
        auto duration_simplekernel = std::chrono::duration_cast<std::chrono::microseconds>(stop_simplekernel - start_simplekernel);
        std::cout << "Simple Kernel gradient time: ";
        std::cout << duration_simplekernel.count();
        std::cout << " microseconds." << std::endl;


        auto start_nofunction = std::chrono::high_resolution_clock::now();

        for(int i = 0 ; i < num_ite; i++)
        {
            double h = std::abs(dis_k(gen)) * 5;
            Eigen::Vector3d xi = Eigen::Vector3d(dis_k(gen), dis_k(gen), dis_k(gen));
            Eigen::Vector3d xj = Eigen::Vector3d(dis_k(gen), dis_k(gen), dis_k(gen));
            learnSPH::kernel::kernel_gradient_with_cubic(xi - xj,h);
            
        }

        auto stop_nofunction = std::chrono::high_resolution_clock::now();
        auto duration_nofunction = std::chrono::duration_cast<std::chrono::microseconds>(stop_nofunction - start_nofunction);
        std::cout << "No function call kernel gradient time: ";
        std::cout << duration_nofunction.count();
        std::cout << " microseconds." << std::endl;


        auto start_nobranch = std::chrono::high_resolution_clock::now();

        for(int i = 0 ; i < num_ite; i++)
        {
            double h = std::abs(dis_k(gen)) * 5;
            Eigen::Vector3d xi = Eigen::Vector3d(dis_k(gen), dis_k(gen), dis_k(gen));
            Eigen::Vector3d xj = Eigen::Vector3d(dis_k(gen), dis_k(gen), dis_k(gen));
            learnSPH::kernel::kernel_gradient_with_cubic_nobranching(xi - xj,h);
            
        }

        auto stop_nobranch = std::chrono::high_resolution_clock::now();
        auto duration_nobranch = std::chrono::duration_cast<std::chrono::microseconds>(stop_nobranch - start_nobranch);
        std::cout << "No function call and no branching kernel gradient time: ";
        std::cout << duration_nobranch.count();
        std::cout << " microseconds." << std::endl;

        return true;
    }
};


TEST_CASE( "Tests for performance", "[performance]" )
{
    KernelFunction kernel;
    CubicSpline cubic;
    GradientCubicSpline gradient;

    //Generate random number for h, and vectors xi, xj
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1, 1);

    //Error variable 
    const double EPSILON = 0.000001;
    int num_ite = 100000;

    SECTION("Testing improvement cubic spline") {
        REQUIRE(cubic.improvement(num_ite, gen, EPSILON));
    }

    SECTION("Testing improvement kernel function") {
        REQUIRE(kernel.improvement(num_ite, gen, EPSILON));

        for(int i = 0 ; i < 1000; i++)
        {
            double h = std::abs(dis(gen)) * 5;
            Eigen::Vector3d xi = Eigen::Vector3d(dis(gen), dis(gen), dis(gen));
            Eigen::Vector3d xj = Eigen::Vector3d(dis(gen), dis(gen), dis(gen));

            REQUIRE(kernel.compare_results(xi,xj,h,EPSILON));
        }
    }

    SECTION("Testing improvement kernel gradient function") {
        REQUIRE(gradient.improvement(num_ite, gen, EPSILON));
    }
}



