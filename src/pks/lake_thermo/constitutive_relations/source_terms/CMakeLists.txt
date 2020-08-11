# -*- mode: cmake -*-

#
#  ATS
#    lake model source terms
#

add_library(lake_thermo_relations_source_terms
  lake_source_terms_evaluator.cc
             )

install(TARGETS lake_thermo_relations_source_terms DESTINATION lib)
