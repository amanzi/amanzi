# -*- mode: cmake -*-

#
#  Amanzi
#   Flow process kernel
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   FLOW_SOURCE_DIR
#   FLOW_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(FLOW)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_SIMPLE_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${NOX_INCLUDE_DIRS})

# Need to move this to roo level -- lpritch
add_definitions("-DUSE_MPI")

#
# Library: flow
#
#add_library(flow MimeticHexLocal.cpp MimeticHex.cpp
#  DiffusionMatrix.cpp DiffusionPrecon.cpp FlowBC.cpp
#  Darcy_PK.cpp DarcyProblem.cpp DarcyMatvec.cpp
#  Richards_PK.cpp RichardsProblem.cpp RichardsNoxInterface.cpp
#  NKA.cpp NKADirFactory.cpp NKADirection.cpp
#  RichardsModelEvaluator.cpp)
file(GLOB flow_source_files "*.cpp")
add_library(flow ${flow_source_files})
add_install_library(flow)
set(flow_tpl_libs ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${NOX_LIBRARIES})
list(REMOVE_DUPLICATES flow_tpl_libs)
target_link_libraries(flow ${flow_tpl_libs})

#
# Install Targets
#
file(GLOB flow_inc_files "*.hpp")
add_install_include_file(${flow_inc_files})

if (BUILD_TESTS) 

    # Add UnitTest include directoy
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Add the STK mesh source directory
    include_directories(${STK_SOURCE_DIR})

    # Copy test directory files if an out of source build
    if (NOT (${FLOW_SOURCE_DIR} EQUAL ${FLOW_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${FLOW_SOURCE_DIR}/test ${FLOW_BINARY_DIR}/test) 
    endif()

    # Add the flow directory to the include paths
    include_directories(${FLOW_SOURCE_DIR})
    include_directories(${SOLVERS_SOURCE_DIR})

    # Test: flow
        
    add_executable(test_flow test/Main.cc test/test_FlowBC.cpp test/test_1d_flow.cpp  test/test_vanGenuchten.cpp)
    target_link_libraries(test_flow flow stk_mesh exodus simple_mesh cell_geometry ${UnitTest_LIBRARIES})
    add_amanzi_test(flow test_flow KIND int)
    add_amanzi_test(flow_parallel_2pe test_flow NPROCS 2 KIND unit)
    add_amanzi_test(flow_parallel_4pe test_flow NPROCS 4 KIND unit)
    add_amanzi_test(flow_parallel_8pe test_flow NPROCS 8 KIND unit)
 

endif()

option(ENABLE_FLOW_EXAMPLES "Build flow examples" OFF)
if (ENABLE_FLOW_EXAMPLES)
   add_subdirectory(examples)
endif()
