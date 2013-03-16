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
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${NOX_INCLUDE_DIRS})
include_directories(${HDF5_C_INCLUDE_DIR})
include_directories(${ASCEMIO_INCLUDE_DIR})

#
# Library: flow
#
set(flow_src_files Flow_PK.cc Flow_State.cc Flow_IO.cc
                   Flow_BC_Factory.cc Flow_BC_Actions.cc Flow_SourceFactory.cc Flow_VandV.cc 
                   Darcy_PK.cc Darcy_TI.cc Darcy_IO.cc Darcy_AztecOO.cc
                   Richards_PK.cc Richards_IO.cc Richards_TI.cc Richards_WRM.cc
                   Richards_SteadyState.cc Richards_Picard.cc Richards_AztecOO.cc
                   Richards_Bundles.cc Richards_Devel.cc Richards_AndersonMixing.cc
                   Matrix_MFD.cc Matrix_MFD_TPFA.cc Matrix_MFD_PLambda.cc Matrix_VAG.cc 
                   WRM_vanGenuchten.cc WRM_BrooksCorey.cc WRM_fake.cc)
set(flow_tpl_libs ${Teuchos_LIBRARIES} ${NOX_LIBRARIES} ${Epetra_LIBRARIES})
if ( Ifpack_ENABLE_HYPRE )
  list(APPEND flow_tpl_libs ${Ifpack_LIBRARIES})
endif()
list(REVERSE flow_tpl_libs)
list(REMOVE_DUPLICATES flow_tpl_libs)
list(REVERSE flow_tpl_libs)
add_amanzi_library(flow SOURCE ${flow_src_files}
                   LINK_LIBS functions mesh_functions ${flow_tpl_libs} time_integration data_structures state)

#
# Library: matrix audit
#
add_amanzi_library(matrix_audit SOURCE Matrix_Audit.cc LINK_LIBS ${flow_tpl_libs})

#
# Install Targets
#
file(GLOB flow_inc_files "*.hh")
add_install_include_file(${flow_inc_files})

if (BUILD_TESTS) 
#if (0)
    # Add UnitTest include directoy
    include_directories(${UnitTest_INCLUDE_DIRS})

    include_directories(${MESH_FACTORY_SOURCE_DIR})

    # Copy test directory files if an out of source build
    if (NOT (${FLOW_SOURCE_DIR} EQUAL ${FLOW_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${FLOW_SOURCE_DIR}/test ${FLOW_BINARY_DIR}/test) 
    endif()

    # Add the flow directory to the include paths
    include_directories(${FLOW_SOURCE_DIR})
    include_directories(${SOLVERS_SOURCE_DIR})

    if ( ENABLE_MSTK_Mesh )
        set(amanzi_libs geometry mesh mesh_audit mesh_factory state flow whetstone data_structures)
 
        # Test: miscaleneous Darcy flow routines
        add_executable(flow_darcy_misc test/Main.cc test/flow_darcy_misc.cc)
        target_link_libraries(flow_darcy_misc ${amanzi_libs} ${UnitTest_LIBRARIES})
        add_amanzi_test(flow_darcy_misc flow_darcy_misc KIND int)
        add_amanzi_test(flow_darcy_parallel_2pe flow_darcy_misc NPROCS 2 KIND unit)
        add_amanzi_test(flow_darcy_parallel_4pe flow_darcy_misc NPROCS 4 KIND unit)

        # Test: steady-state Darcy full flow test
        add_amanzi_test(flow_darcy_source flow_darcy_source
	                KIND int
                        SOURCE test/Main.cc test/flow_darcy_source.cc
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
 
        # Test: transient Darcy flow 
        add_amanzi_test(flow_darcy_2D flow_darcy_2D
	                KIND int
                        SOURCE test/Main.cc test/flow_darcy_2D.cc
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
 
        # Test: transient Darcy flow 
        add_amanzi_test(flow_darcy_dual2D flow_darcy_dual2D
	                KIND int
                        SOURCE test/Main.cc test/flow_darcy_dual2D.cc
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
 
        # Test: transient Darcy flow with a source
        add_amanzi_test(flow_darcy_well flow_darcy_well
	                KIND int
                        SOURCE test/Main.cc test/flow_darcy_well.cc
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

        # Test: miscaleneous Richards flow routines
        add_amanzi_test(flow_richards_convergence flow_richards_convergence
	                KIND int
                        SOURCE test/Main.cc test/flow_richards_convergence.cc
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

        # Test: 2D Richards
        add_amanzi_test(flow_richards_2D flow_richards_2D
	                KIND unit
                        SOURCE test/Main.cc test/flow_richards_2D.cc
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

        add_amanzi_test(flow_richards_seepage flow_richards_seepage
	                KIND unit
                        SOURCE test/Main.cc test/flow_richards_seepage.cc
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

        # Test: pseudo-1D Richards
        add_amanzi_test(flow_richards_bc_cribs flow_richards_bc_cribs
	                KIND unit
                        SOURCE test/Main.cc test/flow_richards_bc_cribs.cc
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
    endif()    
 
    set(amanzi_libs geometry mesh mesh_factory state flow whetstone data_structures)
 
    # Test: van Genuchten flow routines
    add_amanzi_test(flow_vanGenuchten flow_vanGenuchten
                    KIND int
                    SOURCE test/Main.cc test/flow_vanGenuchten.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # Test: flow-bc-factory
    add_amanzi_test(flow_bc_factory flow_bc_factory
                    KIND unit
                    SOURCE test/flow_bc_factory.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
 
    # Test: 3D Richards
    set(amanzi_libs geometry mesh mesh_factory flow state whetstone solvers output data_structures)
    
    add_amanzi_test(flow_richards_tensor flow_richards_tensor
                    KIND unit
                    SOURCE test/Main.cc test/flow_richards_tensor.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    if ( ENABLE_STK_Mesh )
        add_amanzi_test(flow_richards_3D flow_richards_3D
	                KIND unit
                        SOURCE test/Main.cc test/flow_richards_3D.cc
                        LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    endif()
endif()

option(ENABLE_FLOW_EXAMPLES "Build flow examples" OFF)
if (ENABLE_FLOW_EXAMPLES)
   add_subdirectory(examples)
endif()
