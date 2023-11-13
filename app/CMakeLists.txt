cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)

# Local
set(SOURCE_FILES
    main.cpp)

set(DAM_BREAK
    dam_break.cpp)
    
add_executable(learnSPH_app ${SOURCE_FILES})
add_executable(dam_break ${DAM_BREAK})
add_dependencies(learnSPH_app learnSPH)
add_dependencies(dam_break learnSPH)
target_link_libraries(learnSPH_app learnSPH)
target_link_libraries(dam_break learnSPH)
