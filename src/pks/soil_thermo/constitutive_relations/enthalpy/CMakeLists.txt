# -*- mode: cmake -*-

#
#  ATS
#    enthalpy constitutive relations
#

add_library(soil_thermo_relations_enthalpy
  enthalpy_evaluator.cc
             )

install(TARGETS soil_thermo_relations_enthalpy DESTINATION lib)
