# -*- mode: cmake -*-
include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${PROJECT_SOURCE_DIR}/error_handling)
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${MSTK_INCLUDES})
include_directories(${MPI_INCLUDES})
link_directories(${Trilinos_LIBRARY_DIRS})

add_definitions("-DUSE_MPI -DMSTK_HAVE_MPI")

add_library(mstk_mesh Cell_topology.cc Element_category.cc Mesh_maps_mstk.cc)
target_link_libraries(mstk_mesh mesh_data error_handling ${MSTK_LIBS_ALL} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

if (BUILD_TESTS)
  add_executable(test_mstk_mesh 
	test/Main.cc 
	test/test_hex_2x2x2.cc
	test/test_hex_3x3x2.cc
	test/test_hex_4x4x4.cc)
  add_test(mstk_mesh test_mstk_mesh )
  target_link_libraries(test_mstk_mesh mstk_mesh unittest++)

  set(MPI_EXEC_MAX_NUMPROCS 4)
  add_executable(test_parallel_mstk_mesh
        test/Main.cc
        test/test_hex_4x4x4_4P.cc)
  add_test(parallel_mstk_mesh ${MPI_EXEC} ${MPI_EXEC_NUMPROCS_FLAG} ${MPI_EXEC_MAX_NUMPROCS} test_parallel_mstk_mesh)
  target_link_libraries(test_parallel_mstk_mesh mstk_mesh unittest++)

  if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
    execute_process(COMMAND ${CMAKE_COMMAND} -E 
      copy_directory ${CMAKE_SOURCE_DIR}/mstk_mesh/test ${CMAKE_BINARY_DIR}/mstk_mesh/test) 
  endif()
endif()

