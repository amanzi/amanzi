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
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${NOX_INCLUDE_DIRS})
include_directories(${HDF5_C_INCLUDE_DIR})
include_directories(${ASCEMIO_INCLUDE_DIR})
if (ENABLE_CGNS)
    include_directories(${CGNS_INCLUDE_DIRS})
endif()

# Need to move this to roo level -- lpritch
add_definitions("-DUSE_MPI" "-DMSTK_HAVE_MPI")

#
# Library: flow
#
add_library(flow Flow_PK.cpp Flow_State.cpp 
                 Flow_BC_Factory.cpp Flow_Source_Factory.cpp Flow_VandV.cpp 
                 Darcy_PK.cpp Darcy_IO.cpp 
                 Richards_PK.cpp Richards_IO.cpp Richards_TI.cpp Richards_WRM.cpp
                 Richards_SteadyState.cpp Richards_Devel.cpp
                 Matrix_MFD.cpp Matrix_VAG.cpp 
                 WRM_vanGenuchten.cpp WRM_fake.cpp)
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

    set(amanzi_libs geometry mesh mesh_audit mesh_factory simple_mesh mstk_mesh state flow whetstone)
 
    if ( ENABLE_MSTK_Mesh )
        # Test: miscaleneous Darcy flow routines
        add_executable(flow_darcy_misc test/Main.cc test/flow_darcy_misc.cpp)
        target_link_libraries(flow_darcy_misc ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(flow_darcy_misc flow_darcy_misc KIND int)
        add_amanzi_test(flow_darcy_parallel_2pe flow_darcy_misc NPROCS 2 KIND unit)
        add_amanzi_test(flow_darcy_parallel_4pe flow_darcy_misc NPROCS 4 KIND unit)
 
        # Test: transient Darcy flow 
        add_executable(flow_darcy_2D test/Main.cc test/flow_darcy_2D.cpp)
        target_link_libraries(flow_darcy_2D ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(flow_darcy_2D flow_darcy_2D KIND int)
 
        # Test: transient Darcy flow with a source
        add_executable(flow_darcy_well test/Main.cc test/flow_darcy_well.cpp)
        target_link_libraries(flow_darcy_well ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(flow_darcy_well flow_darcy_well KIND int)

        # Test: miscaleneous Richards flow routines
        add_executable(flow_richards_convergence test/Main.cc test/flow_richards_convergence.cpp)
        target_link_libraries(flow_richards_convergence ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(flow_richards_convergence flow_richards_convergence KIND int)
    endif()    
 
    # Test: van Genuchten flow routines
    add_executable(flow_vanGenuchten test/Main.cc test/flow_vanGenuchten.cpp)
    target_link_libraries(flow_vanGenuchten ${amanzi_libs} ${UnitTest_LIBRARIES})
    add_amanzi_test(flow_vanGenuchten flow_vanGenuchten KIND int)

    # Test: flow-bc-factory
    add_executable(flow_bc_factory test/flow_bc_factory.cc)
    target_link_libraries(flow_bc_factory ${amanzi_libs} ${UnitTest_LIBRARIES})
    add_amanzi_test(flow_bc_factory flow_bc_factory KIND unit)
 
    # Test: 2D Richards
    if ( ENABLE_MSTK_Mesh )
        include_directories(${MSTK_INCLUDE_DIRS})
        include_directories(${MSTK_SOURCE_DIR})

        set(amanzi_libs geometry mesh mesh_audit mstk_mesh flow state whetstone solvers output)

        add_executable(flow_richards_2D test/Main.cc test/flow_richards_2D.cpp)
        target_link_libraries(flow_richards_2D ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(flow_richards_2D flow_richards_2D KIND unit)

        add_executable(flow_richards_seepage test/Main.cc test/flow_richards_seepage.cpp)
        target_link_libraries(flow_richards_seepage ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(flow_richards_seepage flow_richards_seepage KIND unit)
    endif()

    # Test: 3D Richards
    set(amanzi_libs geometry mesh mesh_factory mstk_mesh simple_mesh flow state whetstone solvers output)

    add_executable(flow_richards_3D test/Main.cc test/flow_richards_3D.cpp)
    target_link_libraries(flow_richards_3D ${amanzi_libs} ${UnitTest_LIBRARIES})
    add_amanzi_test(flow_richards_3D flow_richards_3D KIND unit)
    
    add_executable(flow_richards_bc_cribs test/Main.cc test/flow_richards_bc_cribs.cpp)
    target_link_libraries(flow_richards_bc_cribs ${amanzi_libs} ${UnitTest_LIBRARIES})
    add_amanzi_test(flow_richards_bc_cribs flow_richards_bc_cribs KIND unit)
    
    add_executable(flow_richards_tensor test/Main.cc test/flow_richards_tensor.cpp)
    target_link_libraries(flow_richards_tensor ${amanzi_libs} ${UnitTest_LIBRARIES})
    add_amanzi_test(flow_richards_tensor flow_richards_tensor KIND unit)
endif()

option(ENABLE_FLOW_EXAMPLES "Build flow examples" OFF)
if (ENABLE_FLOW_EXAMPLES)
   add_subdirectory(examples)
endif()
