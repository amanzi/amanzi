# -*- mode: cmake -*-

#
#  ATS
#    density constitutive relations
#

add_library(soil_thermo_relations_density
  density_evaluator.cc
             )

install(TARGETS soil_thermo_relations_density DESTINATION lib)
