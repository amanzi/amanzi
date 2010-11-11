include_directories(${MPI_INCLUDES})
include_directories(${MOAB_INCLUDES})
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/dbc)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)
include_directories(${PROJECT_SOURCE_DIR}/moab_mesh)
include_directories(${PROJECT_SOURCE_DIR}/moab_mesh)
include_directories(${PROJECT_SOURCE_DIR}/flow)
include_directories(${PROJECT_SOURCE_DIR}/transport)
include_directories(${PROJECT_SOURCE_DIR}/mpc)
include_directories(${PROJECT_SOURCE_DIR}/transport)

# to be elliminated
include_directories(${PROJECT_SOURCE_DIR}/chemistry/includes)
include_directories(${PROJECT_SOURCE_DIR}/flow)

add_definitions("-DUSE_MPI")

add_library(transport Transport_PK.cpp Transport_State.cpp)
target_link_libraries(transport simple_mesh moab_mesh mesh_data dbc mpc flow chemistry ${Trilinos_LIBRARIES})
link_directories(${Trilinos_LIBRARY_DIRS})

if (BUILD_TESTS)
  set(AMANZI_LIBRARIES  simple_mesh moab_mesh mpc flow transport)
  add_executable(test_serial   test/Main.cc test/test_transport.cpp test/test_moab.cpp)
  add_executable(test_parallel test/Main.cc test/test_moab_parallel.cpp)
  add_test(transport test_serial test_parallel)
  target_link_libraries(test_serial   ${AMANZI_LIBRARIES} unittest++ ${Trilinos_LIBRARIES} ${HDF_ALL_LIBS})
  target_link_libraries(test_parallel ${AMANZI_LIBRARIES} unittest++ ${Trilinos_LIBRARIES} ${HDF_ALL_LIBS})
endif()


