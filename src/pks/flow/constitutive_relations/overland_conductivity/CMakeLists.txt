# -*- mode: cmake -*-

#
#  ATS
#    WRM for flow
#

add_library(flow_relations_overland_conductivity
  manning_conductivity_model.cc
  overland_conductivity_evaluator.cc
  unfrozen_fraction_model.cc
  unfrozen_fraction_evaluator.cc
  unfrozen_effective_depth_evaluator.cc
  )

