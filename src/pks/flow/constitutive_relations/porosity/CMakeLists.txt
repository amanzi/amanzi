# -*- mode: cmake -*-

#
#  ATS
#    WRM for flow
#

add_library(flow_relations_porosity
  compressible_porosity_evaluator.cc
  compressible_porosity_model_partition.cc
  )

install(TARGETS flow_relations_porosity DESTINATION lib)

