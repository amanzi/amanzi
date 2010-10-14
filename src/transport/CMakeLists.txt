include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/dbc)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${PROJECT_SOURCE_DIR}/stk_mesh)
include_directories(${PROJECT_SOURCE_DIR}/mpc)
include_directories(${Trilinos_INCLUDE_DIRS})

add_library(transport Transport_PK.cpp)
target_link_libraries(transport stk_mesh dbc mesh_data ${Trilinos_LIBRARIES})
link_directories(${Trilinos_LIBRARY_DIRS})

if (BUILD_TESTS)
  # add additional tests to the end of the following line... eg. test/test_newtest.cpp
  add_executable(test_transport test/Main.cc test/test_example.cpp)
  add_test(transport test_transport)
  # if a new test relies on additional libraries, list them below... eg. flow
  target_link_libraries(test_transport transport stk_mesh unittest++ ${Trilinos_LIBRARIES})
endif()


