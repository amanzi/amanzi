# -*- mode: cmake -*-
include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/error_handling)
include_directories(${PROJECT_SOURCE_DIR}/mpc)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)

include_directories(${Trilinos_INCLUDE_DIRS})

add_library(simple_mesh Mesh_maps_simple.cc)
target_link_libraries(simple_mesh ${Trilinos_LIBRARIES})
link_directories(${Trilinos_LIBRARY_DIRS})

if (BUILD_TESTS)
  add_executable(test_simple_mesh test/Main.cc test/test_maps.cpp test/test_cell_numbering.cc)
  add_test(simple_mesh test_simple_mesh)
  target_link_libraries(test_simple_mesh simple_mesh mpc flow unittest++ ${MPI_LIBS_ALL})
endif()

