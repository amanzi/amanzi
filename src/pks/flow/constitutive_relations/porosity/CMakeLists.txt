# -*- mode: cmake -*-

#
#  ATS
#    WRM for flow
#

add_library(flow_relations_porosity
  compressible_porosity_evaluator.cc
  )


register_evaluator_with_factory(
  HEADERFILE compressible_porosity_evaluator_reg.hh
  LISTNAME   FLOW_CONSTITUTIVE_RELATIONS_POROSITY_EVALUATORS
  )

generate_evaluators_registration_header(
  HEADERFILE flow_constitutive_relations_porosity_evaluators_registration.hh
  LISTNAME   FLOW_CONSTITUTIVE_RELATIONS_POROSITY_EVALUATORS
  INSTALL    True
  )






