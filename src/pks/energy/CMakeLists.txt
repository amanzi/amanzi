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
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${EOS_SOURCE_DIR})
include_directories(${PKS_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${ASCEMIO_INCLUDE_DIR})

#
# Energy registrations
#
register_evaluator_with_factory(
  HEADERFILE twophase_thermal_conductivity_factory_reg.hh
  LISTNAME   ENERGY_CONST_REL_THERMAL_CONDUCTIVITY_REG
  )
register_evaluator_with_factory(
  HEADERFILE twophase_thermal_conductivity_peterslidard_reg.hh
  LISTNAME   ENERGY_CONST_REL_THERMAL_CONDUCTIVITY_REG
  )
register_evaluator_with_factory(
  HEADERFILE twophase_thermal_conductivity_wetdry_reg.hh
  LISTNAME   ENERGY_CONST_REL_THERMAL_CONDUCTIVITY_REG
  )
generate_evaluators_registration_header(
  HEADERFILE energy_thermal_conductivity_registration.hh
  LISTNAME   ENERGY_CONST_REL_THERMAL_CONDUCTIVITY_REG
  INSTALL    True
  )

register_evaluator_with_factory(
  HEADERFILE iem_evaluator_reg.hh
  LISTNAME   ENERGY_INTERNAL_ENERGY_REG
  )
register_evaluator_with_factory(
  HEADERFILE iem_factory_reg.hh
  LISTNAME   ENERGY_INTERNAL_ENERGY_REG
  )
register_evaluator_with_factory(
  HEADERFILE iem_linear_reg.hh
  LISTNAME   ENERGY_INTERNAL_ENERGY_REG
  )
register_evaluator_with_factory(
  HEADERFILE iem_water_vapor_evaluator_reg.hh
  LISTNAME   ENERGY_INTERNAL_ENERGY_REG
  )
generate_evaluators_registration_header(
  HEADERFILE energy_internal_energy_registration.hh
  LISTNAME   ENERGY_INTERNAL_ENERGY_REG
  INSTALL    True
  )

register_evaluator_with_factory(
  HEADERFILE EnergyTwoPhase_PK_reg.hh
  LISTNAME   ENERGY_TWO_PHASE_REG
  )
generate_evaluators_registration_header(
  HEADERFILE energy_twophase_registration.hh
  LISTNAME   ENERGY_TWO_PHASE_REG
  INSTALL    True
  )

#
# Library: energy
#
set(energy_src_files Energy_PK.cc Energy_TI.cc
                     EnergyBoundaryFunction.cc
                     EnergyTwoPhase_PK.cc 
                     Energy_VandV.cc
                     iem_linear.cc
                     iem_water_vapor.cc 
                     iem_water_vapor_evaluator.cc
                     iem_factory.cc iem_evaluator.cc 
                     enthalpy_evaluator.cc
                     twophase_thermal_conductivity_wetdry.cc
                     twophase_thermal_conductivity_peterslidard.cc
                     twophase_energy_evaluator.cc 
                     twophase_thermal_conductivity_evaluator.cc
                     twophase_thermal_conductivity_factory.cc)

set(energy_tpl_libs ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES})
if ( Ifpack_ENABLE_HYPRE )
  list(APPEND energy_tpl_libs ${Ifpack_LIBRARIES})
endif()
add_amanzi_library(energy SOURCE ${energy_src_files}
                   LINK_LIBS functions mesh_functions ${energy_tpl_libs} time_integration data_structures)

#
# Install Targets
#
set(flow_inc_files 
  Energy_PK.hh
  EnergyBoundaryFunction.hh
  EnergyTwoPhase_PK.hh
  energy_twophase_registration.hh
  iem.hh
  iem_linear.hh
  iem_water_vapor.hh
  iem_water_vapor_evaluator.hh
  iem_factory.hh
  iem_evaluator.hh
  energy_internal_energy_registration.hh
  twophase_energy_evaluator.hh
  twophase_thermal_conductivity.hh
  twophase_thermal_conductivity_evaluator.hh
  twophase_thermal_conductivity_factory.hh
  twophase_thermal_conductivity_registration.hh
  twophase_thermal_conductivity_peterslidard.hh
  enthalpy_evaluator.hh)

add_install_include_file(${energy_inc_files})

if (BUILD_TESTS) 

    # Add UnitTest include directoy
    include_directories(${UnitTest_INCLUDE_DIRS})
    include_directories(${MESH_FACTORY_SOURCE_DIR})

    # Copy test directory files if an out of source build
    if (NOT (${ENERGY_SOURCE_DIR} EQUAL ${ENERGY_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${ENERGY_SOURCE_DIR}/test ${ENERGY_BINARY_DIR}/test) 
    endif()

    # Add the flow directory to the include paths
    include_directories(${ENERGY_SOURCE_DIR})

    set(amanzi_libs geometry mesh mesh_factory eos energy eos state operators whetstone data_structures)

    # Test: 1
    add_amanzi_test(energy_matrix_2D energy_matrix_2D
                    KIND int
                    SOURCE test/Main.cc test/energy_matrix_2D.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})

    # Test: 2
    add_amanzi_test(energy_newton_picard energy_newton_picard
                    KIND int
                    SOURCE test/Main.cc test/energy_newton_picard.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
endif()

