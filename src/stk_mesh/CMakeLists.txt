# -*- mode: cmake -*-
include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/error_handling)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${PROJECT_SOURCE_DIR}/exodus)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)

include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})

add_library(stk_mesh Cell_topology.cc Entity_map.cc Mesh_factory.cc
            Mesh.cc Mesh_maps_stk.cc Mesh_maps_stk_hex.cc
            Mesh_maps_stk_exodus.cc Mesh_common.cc)
target_link_libraries(stk_mesh error_handling mesh_data
                      ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
link_directories(${Trilinos_LIBRARY_DIRS})

if (BUILD_TESTS)
  add_executable(test_stk_mesh test/Main.cc test/test_Mesh.cc
                 test/test_Maps.cc test/test_Hex.cc)
  add_test(stk_mesh test_stk_mesh)
  target_link_libraries(test_stk_mesh stk_mesh utils unittest++)

  add_executable(test_stk_mesh_parallel test/Main.cc test/test_Hex.cc)
  target_link_libraries(test_stk_mesh_parallel stk_mesh utils unittest++)
  add_test(stk_mesh_parallel
           ${MPI_EXEC} ${MPI_EXEC_NUMPROCS_FLAG} ${MPI_EXEC_MAX_NUMPROCS} test_stk_mesh_parallel)

  add_executable(test_stk_mesh_read test/Main.cc test/test_Read.cc)
  target_link_libraries(test_stk_mesh_read exodus stk_mesh utils unittest++)
  add_test(stk_mesh_read test_stk_mesh_read)
  add_test(stk_mesh_read_parallel
           ${MPI_EXEC} ${MPI_EXEC_NUMPROCS_FLAG} ${MPI_EXEC_MAX_NUMPROCS} test_stk_mesh_read)

endif()

