include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/dbc)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)
include_directories(${PROJECT_SOURCE_DIR}/chemistry)
include_directories(${PROJECT_SOURCE_DIR}/mpc)
include_directories(${PROJECT_SOURCE_DIR}/gmv)
include_directories(${Trilinos_INCLUDE_DIRS})

add_library(mpc MPC.cpp  State.cpp)
target_link_libraries(mpc simple_mesh gmv mesh_data ${Trilinos_LIBRARIES})
link_directories(${Trilinos_LIBRARY_DIRS})

if (BUILD_TESTS)
  add_executable(test_mpc test/Main.cc test/test_gmv.cpp test/test_driver.cpp)
  add_test(mpc test_mpc)
  target_link_libraries(test_mpc mpc simple_mesh gmv unittest++)
endif()

