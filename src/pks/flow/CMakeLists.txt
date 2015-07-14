# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

add_subdirectory(constitutive_relations)
add_subdirectory(richards)
add_subdirectory(permafrost)
add_subdirectory(overland)
add_subdirectory(icy_overland)
add_subdirectory(overland_head)
add_subdirectory(snow_distribution)

#=========================================================================================

# icy_overland/
register_evaluator_with_factory(
  HEADERFILE icy_overland/icy_height_evaluator_reg.hh
  LISTNAME   FLOW_ICY_OVERLAND_EVAL_FAC
  )

register_evaluator_with_factory(
  HEADERFILE icy_overland/icy_overland_reg.hh
  LISTNAME   FLOW_ICY_OVERLAND_EVAL_FAC
  )

generate_evaluators_registration_header(
  HEADERFILE flow_icy_overland_registration.hh
  LISTNAME   FLOW_ICY_OVERLAND_EVAL_FAC
  INSTALL    True
  )

# overland/
register_evaluator_with_factory(
  HEADERFILE overland/overland_pk_reg.hh
  LISTNAME   FLOW_OVERLAND_HEAD_EVAL_FAC
  )

register_evaluator_with_factory(
  HEADERFILE snow_distribution/snow_distribution_reg.hh
  LISTNAME   FLOW_OVERLAND_HEAD_EVAL_FAC
  )


# overland_head/
register_evaluator_with_factory(
  HEADERFILE overland_head/height_evaluator_reg.hh
  LISTNAME   FLOW_OVERLAND_HEAD_EVAL_FAC
  )

register_evaluator_with_factory(
  HEADERFILE overland_head/effective_height_evaluator_reg.hh
  LISTNAME   FLOW_OVERLAND_HEAD_EVAL_FAC
  )

register_evaluator_with_factory(
  HEADERFILE overland_head/overland_head_pk_reg.hh
  LISTNAME   FLOW_OVERLAND_HEAD_EVAL_FAC
  )

register_evaluator_with_factory(
  HEADERFILE overland_head/overland_head_water_content_evaluator_reg.hh
  LISTNAME   FLOW_OVERLAND_HEAD_EVAL_FAC
  )

generate_evaluators_registration_header(
  HEADERFILE flow_overland_head_registration.hh
  LISTNAME   FLOW_OVERLAND_HEAD_EVAL_FAC
  INSTALL    True
  )


# permafrost/

register_evaluator_with_factory(
  HEADERFILE permafrost/permafrost_pk_reg.hh
  LISTNAME   PERMAFROST_EVAL_FAC
  )

generate_evaluators_registration_header(
  HEADERFILE flow_permafrost_registration.hh
  LISTNAME   PERMAFROST_EVAL_FAC
  INSTALL    True
  )


# richards/

register_evaluator_with_factory(
  HEADERFILE richards/richards_pk_reg.hh
  LISTNAME   RICHARDS_EVAL_FAC
  )

register_evaluator_with_factory(
  HEADERFILE richards/richards_steadystate_reg.hh
  LISTNAME   RICHARDS_EVAL_FAC
  )

generate_evaluators_registration_header(
  HEADERFILE flow_richards_registration.hh
  LISTNAME   RICHARDS_EVAL_FAC
  INSTALL    True
  )





