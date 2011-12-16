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
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${NOX_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${ASCEMIO_INCLUDE_DIR})
if (ENABLE_CGNS)
    include_directories(${CGNS_INCLUDE_DIRS})
endif()

# Need to move this to roo level -- lpritch
add_definitions("-DUSE_MPI")

#
# Library: flow
#
add_library(flow Transient_Richards_PK.cpp Darcy_PK.cpp
                 RichardsProblem.cpp RichardsNoxInterface.cpp RichardsModelEvaluator.cpp
                 DarcyProblem.cpp DarcyMatvec.cpp 
                 DiffusionMatrix.cpp DiffusionPrecon.cpp DiffusionMap.cpp 
                 MimeticHexLocal.cpp MimeticHex.cpp
                 flow-bc-factory.cc vanGenuchtenModel.cpp)
add_install_library(flow)
set(flow_tpl_libs ${Teuchos_LIBRARIES} ${NOX_LIBRARIES} ${Epetra_LIBRARIES})
list(REVERSE flow_tpl_libs)
list(REMOVE_DUPLICATES flow_tpl_libs)
list(REVERSE flow_tpl_libs)
target_link_libraries(flow functions mesh_functions ${flow_tpl_libs} time_integration)
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
    include_directories(${MESH_FACTORY_SOURCE_DIR})

    # Copy test directory files if an out of source build
    if (NOT (${FLOW_SOURCE_DIR} EQUAL ${FLOW_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${FLOW_SOURCE_DIR}/test ${FLOW_BINARY_DIR}/test) 
    endif()

    # Add the flow directory to the include paths
    include_directories(${FLOW_SOURCE_DIR})
    include_directories(${SOLVERS_SOURCE_DIR})

    # Test: miscaleneous flow routines
    add_executable(flow_1D_misc test/Main.cc test/flow_1D_misc.cpp test/flow_vanGenuchten.cpp)
    target_link_libraries(flow_1D_misc flow mesh_factory mesh cell_geometry ${UnitTest_LIBRARIES})
    add_amanzi_test(flow_1D_misc flow_1D_misc KIND int)

    # Test: parallel flow
    add_amanzi_test(flow_parallel_2pe flow_1D_misc NPROCS 2 KIND unit)
    add_amanzi_test(flow_parallel_4pe flow_1D_misc NPROCS 4 KIND unit)
    add_amanzi_test(flow_parallel_8pe flow_1D_misc NPROCS 8 KIND unit)
 
    # Test: flow-bc-factory
    add_executable(flow_bc_factory test/flow_bc_factory.cc)
    target_link_libraries(flow_bc_factory flow mesh_factory mesh error_handling ${UnitTest_LIBRARIES})
    add_amanzi_test(flow_bc_factory flow_bc_factory KIND unit)
 
    # Test: flow-richards
    if (ENABLE_CGNS)
        set(amanzi_libs geometry mesh mesh_audit simple_mesh state flow whetstone solvers output hdf5)

        add_executable(flow_richards_bdf2 test/Main.cc test/flow_richards_bdf2.cpp)
        target_link_libraries(flow_richards_bdf2 ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(flow_richards_bdf2 flow_richards_bdf2 KIND unit)
    endif()
    
endif()

option(ENABLE_FLOW_EXAMPLES "Build flow examples" OFF)
if (ENABLE_FLOW_EXAMPLES)
   add_subdirectory(examples)
endif()
