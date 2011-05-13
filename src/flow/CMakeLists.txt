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
include_directories(${MESH_BASE_SOURCE_DIR})
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
set(flow_tpl_libs ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${NOX_LIBRARIES})
list(REMOVE_DUPLICATES flow_tpl_libs)
target_link_libraries(flow ${flow_tpl_libs})

if (BUILD_TESTS AND ENABLE_MOAB_Mesh)

    # Add UnitTest include directoy
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Copy test directory files if an out of source build
    if (NOT (${FLOW_SOURCE_DIR} EQUAL ${FLOW_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${FLOW_SOURCE_DIR}/test ${FLOW_BINARY_DIR}/test) 
    endif()

    # Add the flow directory to the include paths
    include_directories(${FLOW_SOURCE_DIR})
    include_directories(${SOLVERS_SOURCE_DIR})

    # Test: flow
    # Test depends on MOAB mesh framework
    if ( ENABLE_MOAB_Mesh )
        include_directories(${MOAB_SOURCE_DIR})
        include_directories(${MOAB_INCLUDE_DIRS})
        
        add_executable(test_flow test/Main.cc test/test_FlowBC.cpp test/test_1d_flow.cpp  test/test_vanGenuchten.cpp)
        target_link_libraries(test_flow flow moab_mesh simple_mesh cell_geometry ${UnitTest_LIBRARIES})
        add_amanzi_test(flow test_flow KIND int)
      
    endif()    


endif()

option(ENABLE_FLOW_EXAMPLES "Build flow examples" OFF)
if (ENABLE_FLOW_EXAMPLES)
   add_subdirectory(examples)
endif()
