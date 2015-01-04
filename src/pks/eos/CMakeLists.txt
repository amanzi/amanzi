# -*- mode: cmake -*-

#
#  Amanzi
#    Equations of state
#

#
# Define a project name
# After this command the following varaibles are defined
#   EOS_SOURCE_DIR
#   EOS_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(EOS)

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})

#
# EOS registrations
#
register_evaluator_with_factory(
  HEADERFILE eos_evaluator_reg.hh
  LISTNAME   EOS_EVALUATORS
  )
register_evaluator_with_factory(
  HEADERFILE effective_pressure_evaluator_reg.hh
  LISTNAME   EOS_EVALUATORS
  )
register_evaluator_with_factory(
  HEADERFILE eos_water_reg.hh
  LISTNAME   EOS_EVALUATORS
  )
register_evaluator_with_factory(
  HEADERFILE viscosity_water_reg.hh
  LISTNAME   EOS_EVALUATORS
  )
register_evaluator_with_factory(
  HEADERFILE eos_ideal_gas_reg.hh
  LISTNAME   EOS_EVALUATORS
  )
register_evaluator_with_factory(
  HEADERFILE molar_fraction_gas_evaluator_reg.hh
  LISTNAME   EOS_EVALUATORS
  )
register_evaluator_with_factory(
  HEADERFILE isobaric_eos_evaluator_reg.hh
  LISTNAME   EOS_EVALUATORS
  )
register_evaluator_with_factory(
  HEADERFILE eos_vapor_in_gas_reg.hh
  LISTNAME   EOS_EVALUATORS
  )
register_evaluator_with_factory(
  HEADERFILE vapor_pressure_water_reg.hh
  LISTNAME   EOS_EVALUATORS
  )
generate_evaluators_registration_header(
  HEADERFILE eos_registration.hh
  LISTNAME   EOS_EVALUATORS
  INSTALL    True
  )

#
# Library: eos
#
set(eos_src_files eos_factory.cc 
    eos_evaluator.cc
    molar_fraction_gas_evaluator.cc
    eos_ideal_gas.cc
    eos_vapor_in_gas.cc
    eos_water.cc
    viscosity_water.cc
    viscosity_relation_factory.cc
    vapor_pressure_water.cc
    vapor_pressure_relation_factory.cc
    isobaric_eos_evaluator.cc
    effective_pressure_evaluator.cc
   )

add_amanzi_library(eos
                   SOURCE ${eos_src_files} HEADERS ${eos_inc_files}
                   LINK_LIBS atk)

#
# Install Targets
#
set(eos_inc_files 
    eos.hh
    eos_factory.hh
    eos_evaluator.hh
    eos_registration.hh
    molar_fraction_gas_evaluator.hh
    isobaric_eos_evaluator.hh
    effective_pressure_evaluator.cc
    eos_ideal_gas.hh
    eos_constant_molar_mass.hh
    eos_vapor_in_gas.hh
    eos_water.hh
    viscosity_water.hh
    viscosity_relation.hh
    viscosity_relation_factory.hh
    vapor_pressure_water.hh
    vapor_pressure_relation.hh
    vapor_pressure_relation_factory.hh)

add_install_include_file(${eos_inc_files})

if (BUILD_TESTS) 
endif()

