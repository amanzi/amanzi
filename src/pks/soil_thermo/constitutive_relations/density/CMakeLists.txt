# -*- mode: cmake -*-

#
#  ATS
#    soil density constitutive relations
#

add_library(soil_thermo_relations_density
  soil_density_evaluator.cc
             )

install(TARGETS soil_thermo_relations_density DESTINATION lib)
