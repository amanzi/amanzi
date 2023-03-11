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

# Get the includes directory for chemistry
get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${DBC_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})
include_directories(${ENERGY_SOURCE_DIR})
include_directories(${CHEM_INCLUDES_DIR})
include_directories(${GEOCHEM_SOURCE_DIR})
include_directories(${CHEMPK_SOURCE_DIR})
include_directories(${FLOW_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MULTIPHASE_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${PKS_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${SHALLOW_WATER_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${TRANSPORT_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${ASCEMIO_INCLUDE_DIR})

#
# Library: mpc_tree
#
#
# Install Targets
#
if (ENABLE_AmanziPhysicsModule)
  set(mpc_tree_inc_files
      FractureInsertion.hh
      PK_MPCStrong.hh
      PK_MPCWeak.hh
      PK_MPCSequential.hh
      PK_MPCSubcycled.hh
      PK_MPC.hh
      PK_MPCAdditive.hh
      ChemistryMatrixFracture_PK.hh    
      FlowEnergy_PK.hh
      FlowEnergyMatrixFracture_PK.hh
      FlowMatrixFracture_PK.hh
      MultiphaseMatrixFracture_PK.hh
      TransportMatrixFracture_PK.hh    
      TransportMatrixFractureImplicit_PK.hh
      FlowReactiveTransport_PK.hh
      ReactiveTransport_PK.hh
      ReactiveTransportMatrixFracture_PK.hh
      ShallowWaterTransport_PK.hh
      SurfaceSubsurface_PK.hh)

  set(mpc_tree_src_files 
      FractureInsertion.cc
      PK_MPCWeak.cc 
      PK_MPCSequential.cc 
      PK_MPCSubcycled.cc 
      ChemistryMatrixFracture_PK.cc   
      FlowEnergy_PK.cc
      FlowEnergyMatrixFracture_PK.cc
      FlowMatrixFracture_PK.cc
      MultiphaseMatrixFracture_PK.cc
      TransportMatrixFracture_PK.cc        
      TransportMatrixFractureImplicit_PK.cc
      FlowReactiveTransport_PK.cc
      ReactiveTransport_PK.cc
      ReactiveTransportMatrixFracture_PK.cc
      ShallowWaterTransport_PK.cc
      SurfaceSubsurface_PK.cc)

  set(mpc_tree_link_libs pks flow transport multiphase shallow_water time_integration data_structures state)

else()

  set(mpc_tree_inc_files
      PK_MPCStrong.hh
      PK_MPCWeak.hh
      PK_MPCSequential.hh
      PK_MPCSubcycled.hh
      PK_MPC.hh
      PK_MPCAdditive.hh
      TransportMatrixFracture_PK.hh    
      ReactiveTransport_PK.hh)

  set(mpc_tree_src_files 
      PK_MPCWeak.cc 
      PK_MPCSequential.cc 
      PK_MPCSubcycled.cc 
      TransportMatrixFracture_PK.cc        
      ReactiveTransport_PK.cc)

   set(mpc_tree_link_libs pks transport time_integration data_structures state)

endif()

set(mpc_tree_tpl_libs ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES})
add_amanzi_library(mpc_tree SOURCE ${mpc_tree_src_files}
                   LINK_LIBS ${mpc_tree_link_libs} ${mpc_tree_tpl_libs})


add_install_include_file(${mpc_tree_inc_files} pks_mpc_reg.hh)


if (BUILD_TESTS) 
    # Add UnitTest include directoy
    include_directories(${UnitTest_INCLUDE_DIRS})

    # # Copy test directory files if an out of source build
    # if (NOT (${MPC_TREE_SOURCE_DIR} EQUAL ${MPC_TREE_BINARY_DIR}) )
    #     execute_process(COMMAND ${CMAKE_COMMAND} -E 
    #       copy_directory ${MPC_TREE_SOURCE_DIR}/test ${MPC_TREE_BINARY_DIR}/test) 
    # endif()

    # Add the flow directory to the include paths
    include_directories(${MPC_TREE_SOURCE_DIR})

    set(amanzi_libs geometry mesh mesh_audit mesh_factory state energy operators whetstone data_structures)

    # Test: 1
    #if (ENABLE_AmanziPhysicsModule) 
    #  add_amanzi_test(mpc_pks_test mpc_pks_test
    #                  KIND int
    #                  SOURCE test/Main.cc test/mpc_pks_test.cc
    #                  LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
    #endif()
endif()

