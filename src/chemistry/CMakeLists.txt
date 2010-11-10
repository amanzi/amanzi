include_directories(includes)
include_directories(${PROJECT_SOURCE_DIR}/transport)
include_directories(${PROJECT_SOURCE_DIR}/flow)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)
include_directories(${PROJECT_SOURCE_DIR}/mpc)

add_subdirectory(utility)
add_subdirectory(solvers)
add_subdirectory(base_chemistry)
add_subdirectory(reactions)
add_subdirectory(process_kernel)
add_subdirectory(test_drivers)

include_directories(${PROJECT_SOURCE_DIR}/mpc)
include_directories(${PROJECT_SOURCE_DIR}/transport)
include_directories(${PROJECT_SOURCE_DIR}/flow)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)
include_directories(${PROJECT_SOURCE_DIR}/util)
include_directories(${Trilinos_INCLUDE_DIRS})

link_directories(${Trilinos_LIBRARY_DIRS})

if (BUILD_TESTS)
  add_executable(test_chemistry test/Main.cc test/test_driver.cpp)
  add_test(chemistry test_chemistry)
  target_link_libraries(test_chemistry mpc simple_mesh chemistry transport flow utils unittest++)
endif()

