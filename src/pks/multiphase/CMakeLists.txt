# -*- mode: cmake -*-

#
#  Amanzi
#   Multiphase process kernel
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   MULTIPHASE_SOURCE_DIR
#   MULTIPHASE_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MULTIPHASE)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${PLATFORM_INTERFACE_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${PKS_SOURCE_DIR})
include_directories(${MPC_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})
include_directories(${FLOW_SOURCE_DIR})
include_directories(${MULTIPHASE_SOURCE_DIR})


# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${HDF5_C_INCLUDE_DIR})
include_directories(${ASCEMIO_INCLUDE_DIR})

#
# Multiphase registrations
#
register_evaluator_with_factory( HEADERFILE Pressure_PK_reg.hh 
  Saturation_PK_reg.hh
  Phase2_PK_reg.hh 
  Phase1_PK_reg.hh
  Comp_PK_reg.hh
  MPCoupled_PK_reg.hh
  MPMC_PK_reg.hh
  Phase_Constraint_PK_reg.hh
  Reduced2p2c_PK_reg.hh
	LISTNAME   MULTIPHASE_REGISTRATIONS )

generate_evaluators_registration_header( 
  HEADERFILE pks_pressure_registration.hh 
  pks_saturation_registration.hh
  pks_phase2_registration.hh 
  pks_phase1_registration.hh
  pks_comp_registration.hh
  pks_mpmc_registration.hh
  pks_phase_constraint_registration.hh
  pks_mpcoupled_registration.hh
  LISTNAME   MULTIPHASE_REGISTRATIONS
  INSTALL    True )


#
# Library: multiphase 
#
set(multiphase_src_files Pressure_PK.cc Pressure_PK_TI.cc Pressure_PK_IO.cc
                         Saturation_PK.cc Saturation_PK_TI.cc Saturation_PK_IO.cc
                         Phase2_PK.cc Phase2_PK_TI.cc Phase2_PK_IO.cc
                         Phase1_PK.cc Phase1_PK_TI.cc Phase1_PK_IO.cc
                         TotalPhase_PK.cc TotalPhase_PK_TI.cc
                         Comp_PK.cc Comp_PK_TI.cc Comp_PK_IO.cc
                         Phase_Constraint_PK.cc
                         MPMC_PK.cc
                         CompW.cc CompW_TI.cc CompH.cc CompH_TI.cc
                         GasConstraint.cc
                         Reduced2p2c_PK.cc
                         MPCoeff.cc
                         RelativePermeability_IO.cc FractionalFlow.cc
                         RelativePermeability.cc WRM_Simple.cc
                         WRM_VanGenuchten.cc WRM_BrooksCorey.cc
                         Multiphase_BC_Factory.cc CapillaryPressure.cc
                         MPCoupled_PK.cc
                         )
#set(flow_tpl_libs ${Teuchos_LIBRARIES} ${NOX_LIBRARIES} ${Epetra_LIBRARIES})
#list(REVERSE flow_tpl_libs)
#list(REMOVE_DUPLICATES flow_tpl_libs)
#list(REVERSE flow_tpl_libs)
add_amanzi_library(multiphase 
                   SOURCE ${multiphase_src_files}
                   LINK_LIBS)
#
# Install Targets
#
#file(GLOB flow_inc_files "*.hh")
set(multiphase_inc_files 
	Pressure_PK.hh
  Saturation_PK.hh
  Phase2_PK.hh
  Phase1_PK.hh
  TotalPhase_PK.hh
  Comp_PK.hh
  CompW.hh
  CompH.hh
  Phase_Constraint_PK.hh
  GasConstraint.hh
  MPMC_PK.hh
  Reduced2p2c_PK.hh
  MPCoeff.hh
  RelativePermeability.hh
  CapillaryPressure.hh
  FractionalFlow.hh
  WaterRetentionModel.hh
  WRM_Simple.hh
  WRM_VanGenuchten.hh
  WRM_BrooksCorey.hh
  Multiphase_BC_Factory.hh
  MPCoupled_PK.hh
  )  


add_install_include_file(${multiphase_inc_files})

if (BUILD_TESTS) 
  include_directories(${MESH_FACTORY_SOURCE_DIR})

  # Add UnitTest includes
  include_directories(${UnitTest_INCLUDE_DIRS})

  # Add MPC include directory. Need to remove this required path. -- lpritch
  include_directories(${MULTIPHASE_SOURCE_DIR})

    # Copy test subdirectory for out of source builds
	if (NOT (MULTIPHASE_SOURCE_DIR STREQUAL MULTIPHASE_BINARY_DIR))
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
			copy_directory ${MULTIPHASE_SOURCE_DIR}/test ${MULTIPHASE_BINARY_DIR}/test) 
    endif()

  set(amanzi_libs multiphase state mesh_factory whetstone operators flow
                    geometry mesh_functions functions pks)

  # Unit tests:
if(0)
  add_amanzi_test(test_pressure_pk test_pressure_pk
                  KIND      int
                  SOURCE    test/Main.cc test/test_pressure_pk.cc
		  LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

  add_amanzi_test(test_pressure_pk_functional test_pressure_pk_functional
                  KIND      int
                  SOURCE    test/Main.cc test/test_pressure_pk_functional.cc
      LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

  add_amanzi_test(test_pressure_convergence test_pressure_convergence
                  KIND      int
                  SOURCE    test/Main.cc test/test_pressure_convergence.cc
      LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

  add_amanzi_test(test_saturation_bl test_saturation_bl
                  KIND      int
                  SOURCE    test/Main.cc test/test_saturation_bl.cc
      LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

  add_amanzi_test(test_sequential test_sequential
                  KIND      int
                  SOURCE    test/Main.cc test/test_sequential.cc
      LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
      
  add_amanzi_test(test_seq_q5 test_seq_q5
                  KIND      int
                  SOURCE    test/Main.cc test/test_seq_q5.cc
      LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

  add_amanzi_test(test_seq_q5_3d test_seq_q5_3d
                  KIND      int
                  SOURCE    test/Main.cc test/test_seq_q5_3d.cc
      LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

  add_amanzi_test(test_phase1_pk test_phase1_pk
                  KIND      int
                  SOURCE    test/Main.cc test/test_phase1_pk.cc
      LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

  add_amanzi_test(test_phase2_pk test_phase2_pk
                  KIND      int
                  SOURCE    test/Main.cc test/test_phase2_pk.cc
      LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
endif()

  add_amanzi_test(test_mpcoupled_q5 test_mpcoupled_q5
                  KIND      int
                  SOURCE    test/Main.cc test/test_mpcoupled_q5.cc
      LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

if(0)
  add_amanzi_test(test_mpcoupled test_mpcoupled
                  KIND      int
                  SOURCE    test/Main.cc test/test_mpcoupled.cc
      LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

  add_amanzi_test(test_mpcoupled_pc test_mpcoupled_pc
                  KIND      int
                  SOURCE    test/Main.cc test/test_mpcoupled_pc.cc
      LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

  add_amanzi_test(test_mpcoupled_3d test_mpcoupled_3d
                  KIND      int
                  SOURCE    test/Main.cc test/test_mpcoupled_3d.cc
      LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

  add_amanzi_test(test_constraint test_constraint
                  KIND      int
                  SOURCE    test/Main.cc test/test_constraint.cc
      LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

  add_amanzi_test(test_comp_pk test_comp_pk
                  KIND      int
                  SOURCE    test/Main.cc test/test_comp_pk.cc
      LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

  add_amanzi_test(test_mpmc test_mpmc
                  KIND      int
                  SOURCE    test/Main.cc test/test_mpmc.cc
      LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
endif()
  add_amanzi_test(test_reduced_2p2c test_reduced_2p2c
                  KIND      int
                  SOURCE    test/Main.cc test/test_reduced_2p2c.cc
      LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
endif()

