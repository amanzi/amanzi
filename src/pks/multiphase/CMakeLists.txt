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
# Library: multiphase 
#
set(multiphase_src_files Multiphase_PK.cc Multiphase_TI.cc
                         CompW.cc CompW_TI.cc CompH.cc CompH_TI.cc
                         GasConstraint.cc
                         Reduced2p2c_PK.cc
                         MPCoeff.cc
                         WRMmp_Simple.cc WRMmp_VanGenuchten.cc WRMmp_BrooksCorey.cc
                         CapillaryPressure.cc
                         )

add_amanzi_library(multiphase 
                   SOURCE ${multiphase_src_files}
                   LINK_LIBS pks operators state functions 
                             mesh_functions time_integration data_structures
                             flow)

#
# Install Targets
#
set(multiphase_inc_files 
  Multipahase_PK.hh
  ModelMeshPartition.hh
  CompW.hh
  CompH.hh
  GasConstraint.hh
  Reduced2p2c_PK.hh
  MPCoeff.hh
  CapillaryPressure.hh
  WRMmp.hh
  WRMmp_Simple.hh
  WRMmp_VanGenuchten.hh
  WRMmp_BrooksCorey.hh)  

add_install_include_file(${multiphase_inc_files})

if (BUILD_TESTS) 
  include_directories(${MESH_FACTORY_SOURCE_DIR})

  # Add UnitTest includes
  include_directories(${UnitTest_INCLUDE_DIRS})

  # Add MPC include directory. Need to remove this required path. -- lpritch
  include_directories(${MULTIPHASE_SOURCE_DIR})

  # Copy test subdirectory for out of source builds
  if (NOT (MULTIPHASE_SOURCE_DIR STREQUAL MULTIPHASE_BINARY_DIR))
    file(GLOB DataFiles "${MULTIPHASE_SOURCE_DIR}/test/*.xml"
                        "${MULTIPHASE_SOURCE_DIR}/test/*.exo")
    file(COPY ${DataFiles} DESTINATION ${MULTIPHASE_BINARY_DIR}/test/)
  endif()

  set(amanzi_libs multiphase state mesh_factory whetstone operators flow
                  geometry mesh_functions functions pks)

  # Unit tests:
  add_amanzi_test(test_reduced_2p2c test_reduced_2p2c
                  KIND      int
                  SOURCE    test/Main.cc test/test_reduced_2p2c.cc
                  LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

  add_amanzi_test(multiphase_model1 multipase_model1
                  KIND      int
                  SOURCE    test/Main.cc test/multiphase_model1.cc
                  LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
endif()

