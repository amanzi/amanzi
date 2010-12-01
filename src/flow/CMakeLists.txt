# -*- mode: cmake -*-
include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/error_handling)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)
include_directories(${PROJECT_SOURCE_DIR}/moab_mesh)
include_directories(${PROJECT_SOURCE_DIR}/mpc)
include_directories(${PROJECT_SOURCE_DIR}/flow)
include_directories(${PROJECT_SOURCE_DIR}/gmv)
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDES})
add_definitions("-DUSE_MPI")

add_library(flow Flow_PK.cpp MimeticHexLocal.cpp MimeticHex.cpp
  DiffusionMatrix.cpp DiffusionPrecon.cpp DarcyProblem.cpp
  DarcyMatvec.cpp FlowBC.cpp)
# add_library(flow upper_packed_matrix.cpp mimetic_hex.cpp Flow_PK.cpp
#  Darcy_problem.cpp Flow_BCs.cpp cell_geometry.cpp DiffusionMatrix.cpp
#  DiffusionPrecon.cpp DarcyProblem.cpp DarcyMatvec.cpp FlowBC.cpp)
target_link_libraries(flow ${MESH_LIBRARY} ${Trilinos_LIBRARIES})
link_directories(${Trilinos_LIBRARY_DIRS})

if (BUILD_TESTS)
#  add_executable(test_1d_flow test/test_1d_flow.cpp)
#  target_link_libraries(test_1d_flow flow simple_mesh utils)
#  add_test(flow test_1d_flow)

#  add_executable(test_flow test/Main.cc test/test_FlowBC.cpp test/test_1d_flow.cpp)
  add_executable(test_flow test/Main.cc test/test_FlowBC.cpp )
#  add_executable(test_flow test/Main.cc test/test_FlowBC.cpp test_mimetic.cpp)
  add_test(flow test_flow)
  if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
    execute_process(COMMAND ${CMAKE_COMMAND} -E 
      copy_directory ${CMAKE_SOURCE_DIR}/flow/test ${CMAKE_BINARY_DIR}/flow/test) 
  endif()
  target_link_libraries(test_flow flow simple_mesh ${MESH_LIBRARY} utils unittest++ ${Trilinos_LIBRARIES} ${MPI_LIBS_ALL})

endif()

