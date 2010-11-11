include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${PROJECT_SOURCE_DIR}/dbc)
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDES})
include_directories(${MPI_INCLUDES})
link_directories(${Trilinos_LIBRARY_DIRS})

add_definitions("-DUSE_MPI")

add_library(moab_mesh Cell_topology.cc Element_category.cc Mesh_maps_moab.cc)
target_link_libraries(moab_mesh mesh_data dbc ${MOAB_LIBS_ALL} ${Trilinos_LIBRARIES})

if (BUILD_TESTS)
  add_executable(test_moab_mesh 
	test/Main.cc 
	test/test_hex_2x2x2.cc
	test/test_hex_3x3x2.cc
	test/test_hex_4x4x4.cc)
  add_test(moab_mesh test_moab_mesh )
  if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
    execute_process(COMMAND ${CMAKE_COMMAND} -E 
      copy_directory ${CMAKE_SOURCE_DIR}/moab_mesh/test ${CMAKE_BINARY_DIR}/moab_mesh/test) 
  endif()
  target_link_libraries(test_moab_mesh moab_mesh unittest++ ${MOAB_LIBS_ALL})
endif()

