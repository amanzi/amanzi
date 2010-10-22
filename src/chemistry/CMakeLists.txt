project(ascem-geochemistry CXX)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -ansi -Wextra")

cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE debug)

include_directories(./includes)

add_subdirectory(utility)
add_subdirectory(solvers)
add_subdirectory(base_chemistry)
add_subdirectory(reactions)
add_subdirectory(test_drivers)
