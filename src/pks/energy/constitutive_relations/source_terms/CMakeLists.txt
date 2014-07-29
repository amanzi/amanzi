# -*- mode: cmake -*-

#
#  ATS
#    internal energy constitutive relations
#

add_library(energy_relations_source_terms
  advected_energy_source_evaluator.cc
             )

install(TARGETS energy_relations_source_terms DESTINATION lib)
