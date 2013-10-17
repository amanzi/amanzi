# -*- mode: cmake -*-

#
#  ATS
#    WRM for flow
#

add_subdirectory(models)

include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/wrm/models)

add_library(flow_relations_wrm
  wrm_partition.cc
  wrm_permafrost_evaluator.cc
  wrm_evaluator.cc
  rel_perm_evaluator.cc
  pc_liquid_evaluator.cc
  pc_ice_evaluator.cc
  )


register_evaluator_with_factory(
  HEADERFILE pc_ice_evaluator_reg.hh
  LISTNAME   WMR_EVALUATORS_FACTORIES
  )

register_evaluator_with_factory(
  HEADERFILE pc_liquid_evaluator_reg.hh
  LISTNAME   WMR_EVALUATORS_FACTORIES
  )

register_evaluator_with_factory(
  HEADERFILE wrm_evaluator_reg.hh
  LISTNAME   WMR_EVALUATORS_FACTORIES
  )

register_evaluator_with_factory(
  HEADERFILE wrm_permafrost_evaluator_reg.hh
  LISTNAME   WMR_EVALUATORS_FACTORIES
  )

generate_evaluators_registration_header(
  HEADERFILE wrm_evaluators_factories_registration.hh
  LISTNAME   WMR_EVALUATORS_FACTORIES
  INSTALL    True
  )



