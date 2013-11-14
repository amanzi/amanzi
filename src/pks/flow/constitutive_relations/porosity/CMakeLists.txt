# -*- mode: cmake -*-

#
#  ATS
#    WRM for flow
#

add_library(flow_relations_porosity
  compressible_porosity_evaluator.cc
  )

install(TARGETS flow_relations_porosity DESTINATION lib)

