# -*- mode: cmake -*-
include_directories(${PROJECT_SOURCE_DIR}/utils)
include_directories(${PROJECT_SOURCE_DIR}/error_handling)
include_directories(${PROJECT_SOURCE_DIR}/mesh)
include_directories(${PROJECT_SOURCE_DIR}/mesh_data)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)
include_directories(${PROJECT_SOURCE_DIR}/moab_mesh)
include_directories(${PROJECT_SOURCE_DIR}/mpc)
include_directories(${PROJECT_SOURCE_DIR}/flow)
include_directories(${PROJECT_SOURCE_DIR}/gmv)
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDES})
add_definitions("-DUSE_MPI")

add_library(flow MimeticHexLocal.cpp MimeticHex.cpp
  DiffusionMatrix.cpp DiffusionPrecon.cpp FlowBC.cpp
  Darcy_PK.cpp DarcyProblem.cpp DarcyMatvec.cpp
  Richards_PK.cpp RichardsProblem.cpp RichardsNoxInterface.cpp
  NKA.cpp NKADirFactory.cpp NKADirection.cpp
  RichardsModelEvaluator.cpp)
target_link_libraries(flow ${MESH_LIBRARY} ${Trilinos_LIBRARIES})
link_directories(${Trilinos_LIBRARY_DIRS})

if (BUILD_TESTS)
  add_executable(test_flow test/Main.cc test/test_FlowBC.cpp test/test_1d_flow.cpp)
  add_test(flow test_flow)
  target_link_libraries(test_flow flow moab_mesh simple_mesh ${MESH_LIBRARY} cell_geometry unittest++ ${Trilinos_LIBRARIES})
  if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
    execute_process(COMMAND ${CMAKE_COMMAND} -E 
      copy_directory ${CMAKE_SOURCE_DIR}/flow/test ${CMAKE_BINARY_DIR}/flow/test) 
  endif()
endif()

OPTION(ENABLE_FLOW_EXAMPLES "Build flow examples" OFF)
if (ENABLE_FLOW_EXAMPLES)
 add_subdirectory(examples)
endif ()
