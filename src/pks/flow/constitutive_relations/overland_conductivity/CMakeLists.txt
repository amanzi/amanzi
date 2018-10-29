# -*- mode: cmake -*-

#
#  ATS
#    WRM for flow
#

add_library(flow_relations_overland_conductivity
  manning_conductivity_model.cc
  split_denominator_conductivity_model.cc
  ponded_depth_passthrough_conductivity_model.cc
  overland_conductivity_evaluator.cc
  unfrozen_fraction_model.cc
  unfrozen_fraction_evaluator.cc
  unfrozen_effective_depth_evaluator.cc
  surface_relperm_evaluator.cc
  surface_relperm_model_factory.cc
  unfrozen_fraction_relperm_model.cc
  zero_uf_relperm_model.cc
  one_uf_relperm_model.cc
  fractional_conductance_evaluator.cc	
  )

install(TARGETS flow_relations_overland_conductivity DESTINATION lib)

