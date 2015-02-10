# -*- mode: cmake -*-

#
#  Amanzi
#   Energy process kernel
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   MPC_TREE_SOURCE_DIR
#   MPC_TREE_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MPC_TREE)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})
include_directories(${PKS_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${ASCEMIO_INCLUDE_DIR})

#
# PK registrations
#
register_evaluator_with_factory( HEADERFILE FlowReactiveTransport_PK_reg.hh 
  LISTNAME   PKS_FRT_REGISTRATIONS )

generate_evaluators_registration_header( 
  HEADERFILE pks_flowreactivetransport_registration.hh 
  LISTNAME   PKS_FRT_REGISTRATIONS
  INSTALL    True )


#
# Library: energy
#
set(mpc_tree_src_files MPCWeak.cc MPCSubcycled.cc FlowReactiveTransport_PK.cc)
set(mpc_tree_tpl_libs ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES})
add_amanzi_library(mpc_tree SOURCE ${mpc_tree_src_files}
                   LINK_LIBS ${mpc_tree_tpl_libs} time_integration data_structures state)

#
# Install Targets
#
set(mpc_tree_inc_files MPCStrong.hh MPCWeak.hh MPCSubcycled.hh MPC_tmp.hh FlowReactiveTransport_PK.hh)

add_install_include_file(${mpc_tree_inc_files})

if (BUILD_TESTS) 
  
#if (0)
    # Add UnitTest include directoy
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Copy test directory files if an out of source build
    if (NOT (${MPC_TREE_SOURCE_DIR} EQUAL ${MPC_TREE_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${MPC_TREE_SOURCE_DIR}/test ${MPC_TREE_BINARY_DIR}/test) 
    endif()

    # Add the flow directory to the include paths
    include_directories(${MPC_TREE_SOURCE_DIR})

    set(amanzi_libs geometry mesh mesh_audit mesh_factory state energy operators whetstone data_structures)

    # Test: 1
    #add_amanzi_test(energy_test energy_test
    #                KIND int
    #                SOURCE test/Main.cc test/energy_test.cc
    #                LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
endif()

