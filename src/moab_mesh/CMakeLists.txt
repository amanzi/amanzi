include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/dbc)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIR})
include_directories(${MPI_INCLUDE_DIR})
link_directories(${MOAB_LIB_DIR} ${Trilinos_LIBRARY_DIRS})

add_definitions("-DUSE_MPI")

add_library(moab_mesh Cell_topology.cc Element_category.cc Mesh_maps.cc)
target_link_libraries(moab_mesh dbc mesh_data ${MOAB_LIBS_ALL} ${Trilinos_LIBRARIES})

if (BUILD_TESTS)
  add_executable(test_moab_mesh 
	test/Main.cc 
	test/test_hex_2x2x2.cc
	test/test_hex_3x3x2.cc)
  add_test(moab_mesh test_moab_mesh)
  target_link_libraries(test_moab_mesh moab_mesh unittest++)
endif()

