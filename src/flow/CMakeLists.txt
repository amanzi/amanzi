include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/dbc)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${PROJECT_SOURCE_DIR}/stk_mesh)
include_directories(${PROJECT_SOURCE_DIR}/mpc)
include_directories(${Trilinos_INCLUDE_DIRS})

add_library(flow upper_packed_matrix.cpp mimetic_hex.cpp Flow_PK.cpp Darcy_problem.cpp Flow_BCs.cpp cell_geometry.cpp)
target_link_libraries(flow stk_mesh ${Trilinos_LIBRARIES})
link_directories(${Trilinos_LIBRARY_DIRS})

if (BUILD_TESTS)
  add_executable(test_flow test/Main.cc test/test_UPM.cpp test/test_mimetic.cpp test/test_Flow_BCs.cpp)
#  add_executable(test_flow test/Main.cc test/test_UPM.cpp test/test_mimetic.cpp)
  add_test(flow test_flow)
  target_link_libraries(test_flow flow stk_mesh unittest++ ${Trilinos_LIBRARIES})
endif()

