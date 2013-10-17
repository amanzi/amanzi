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

register_evaluator_with_factory(
  HEADERFILE unfrozen_effective_depth_evaluator_reg.hh
  LISTNAME   FLOW_CONSTITUTIVE_RELATIONS_OVERLAND_CONDUCTIVITY
  )

register_evaluator_with_factory(
  HEADERFILE unfrozen_fraction_evaluator_reg.hh
  LISTNAME   FLOW_CONSTITUTIVE_RELATIONS_OVERLAND_CONDUCTIVITY
  )

generate_evaluators_registration_header(
  HEADERFILE flow_constitutive_relations_overland_conductivity_registration.hh
  LISTNAME   FLOW_CONSTITUTIVE_RELATIONS_OVERLAND_CONDUCTIVITY
  INSTALL    True
  )







