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
#   ENERGY_SOURCE_DIR
#   ENERGY_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(ENERGY)

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${DBC_SOURCE_DIR})
include_directories(${EOS_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${PKS_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})

#
# Library: energy
#
set(energy_src_files Energy_PK.cc
                     EnergyOnePhase_PK.cc EnergyOnePhase_TI.cc
                     EnergyTwoPhase_PK.cc EnergyTwoPhase_TI.cc
                     Energy_VandV.cc
                     IEM_Linear.cc
                     IEM_WaterVapor.cc 
                     IEM_WaterVaporEvaluator.cc
                     IEMFactory.cc IEMEvaluator.cc 
                     EnthalpyEvaluator.cc
                     TCM_WetDry_TwoPhase.cc
                     TCM_PetersLidard_TwoPhase.cc
                     TCMEvaluator_OnePhase.cc
                     TCMEvaluator_TwoPhase.cc
                     TCMFactory_TwoPhase.cc
                     TotalEnergyEvaluator.cc)

set(energy_tpl_libs ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES})
if (Ifpack_ENABLE_HYPRE)
  list(APPEND energy_tpl_libs ${Ifpack_LIBRARIES} ${HYPRE_LIBRARY})
endif()
add_amanzi_library(energy SOURCE ${energy_src_files}
                   LINK_LIBS operators eos pks functions mesh_functions
                             time_integration data_structures
                             ${energy_tpl_libs})

#
# Install Targets
#
set(flow_inc_files 
  Energy_PK.hh
  EnergyDefs.hh
  EnergyOnePhase_PK.hh
  EnergyTwoPhase_PK.hh
  pks_energy_registration.hh
  IEM.hh
  IEM_Linear.hh
  IEM_Tabular.hh
  IEM_WaterVapor.hh
  IEM_WaterVaporEvaluator.hh
  IEMFactory.hh
  IEMEvaluator.hh
  energy_iem_registration.hh
  energy_tcm_registration.hh
  TotalEnergyEvaluator.hh
  TCM_PetersLidard_TwoPhase.hh
  TCM_TwoPhase.hh
  TCMEvaluator_OnePhase.hh
  TCMEvaluator_TwoPhase.hh
  TCMFactory_TwoPhase.hh
  EnthalpyEvaluator.hh)

add_install_include_file(${energy_inc_files})

if (BUILD_TESTS) 

    # Add UnitTest include directoy
    include_directories(${UnitTest_INCLUDE_DIRS})
    include_directories(${MESH_FACTORY_SOURCE_DIR})

    # Copy test directory files if an out of source build
    if (NOT (${ENERGY_SOURCE_DIR} EQUAL ${ENERGY_BINARY_DIR}) )
      file(GLOB DataFiles "${ENERGY_SOURCE_DIR}/test/*.xml")
      file(COPY ${DataFiles} DESTINATION ${ENERGY_BINARY_DIR}/test/)
    endif()

    # Add the flow directory to the include paths
    include_directories(${ENERGY_SOURCE_DIR})

    set(amanzi_libs geometry mesh mesh_factory eos energy pks eos state operators whetstone data_structures)

    # Test: 1
    add_amanzi_test(energy_matrix_2D energy_matrix_2D
                    KIND int
                    SOURCE test/Main.cc test/energy_matrix_2D.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: 2
    add_amanzi_test(energy_newton_picard energy_newton_picard
                    KIND int
                    SOURCE test/Main.cc test/energy_newton_picard.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: 3
    add_amanzi_test(energy_convergence energy_convergence
                    KIND int
                    SOURCE test/Main.cc test/energy_convergence.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: 4
    add_amanzi_test(energy_one_phase energy_one_phase
                    KIND int
                    SOURCE test/Main.cc test/energy_one_phase.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # Test: 5
    add_amanzi_test(energy_convergence_src energy_convergence_src
                    KIND int
                    SOURCE test/Main.cc test/energy_convergence_src.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})
                    
endif()

