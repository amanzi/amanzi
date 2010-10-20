include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/dbc)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)
include_directories(${PROJECT_SOURCE_DIR}/mpc)
include_directories(${Trilinos_INCLUDE_DIRS})

add_library(transport Transport_PK.cpp Transport_State.cpp)
target_link_libraries(transport simple_mesh mesh_data dbc mpc ${Trilinos_LIBRARIES})
link_directories(${Trilinos_LIBRARY_DIRS})

if (BUILD_TESTS)
  # add additional tests to the end of the following line... eg. test/test_newtest.cpp
  add_executable(test_transport test/Main.cc test/test_example.cpp)
  add_test(transport test_transport)
  # if a new test relies on additional libraries, list them below... eg. flow
  target_link_libraries(test_transport transport simple_mesh unittest++ ${Trilinos_LIBRARIES})
endif()


