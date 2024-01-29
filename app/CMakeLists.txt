cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)

# Local
set(SOURCE_FILES
    main.cpp)

# set(BENCHMARK_FILES
#     mcubes_benchmark.cpp)
    
add_executable(learnSPH_app ${SOURCE_FILES})
add_dependencies(learnSPH_app learnSPH)
target_link_libraries(learnSPH_app learnSPH)

# add_executable(learnSPH_benchmarkMC ${BENCHMARK_FILES})
# add_dependencies(learnSPH_benchmarkMC learnSPH)
# target_link_libraries(learnSPH_benchmarkMC learnSPH)
