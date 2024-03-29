# External dependencies
## OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

## Eigen
set( EIGEN3_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/extern/Eigen )
if( NOT EIGEN3_INCLUDE_DIR )
    message( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
endif()


set(SOURCE_FILES
	io.h
	io.cpp
	sampling.h
	sampling.cpp
	kernel.cpp
	kernel.h
	neighborhood_search.cpp
	neighborhood_search.h
	time_integration.cpp
	time_integration.h
	types.h
	acceleration.h
	acceleration.cpp
	densities.h
	densities.cpp
	geometry.h
	geometry.cpp
	utils.h
	utils.cpp
	simulations_setup.h
	simulations_setup.cpp
	marching_cubes_lut.h
	marching_cubes.h
	marching_cubes.cpp
	theta_functions.h
	theta_functions.cpp
	pbf.h
	pbf.cpp
	emitter.h
	emitter.cpp
	surface_tension.h
	surface_tension.cpp)

    
# Create a library version so it can be linked with other modules
add_library(learnSPH ${SOURCE_FILES})
target_link_libraries(learnSPH CompactNSearch)
target_include_directories(learnSPH PUBLIC ${EIGEN3_INCLUDE_DIR} )
target_include_directories(learnSPH PUBLIC ${PROJECT_SOURCE_DIR}/extern ) # So you can include header files under learnSPH/extern
# target_link_libraries(learnSPH PRIVATE nlohmann_json::nlohmann_json)