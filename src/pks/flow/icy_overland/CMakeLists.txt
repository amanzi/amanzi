# -*- mode: cmake -*-
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/factory)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/bc_factory)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/overland_conductivity)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/elevation)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad/upwind_scheme)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/overland_head)

add_library(pk_flow_icy_overland
    icy_overland.cc
    icy_height_evaluator.cc
    overland_head_icy_water_content_evaluator.cc)


register_evaluator_with_factory(
  HEADERFILE icy_overland_reg.hh
  LISTNAME   FLOW_ICY_OVERLAND_EVAL_FAC
  )

generate_evaluators_registration_header(
  HEADERFILE icy_overland_evaluators_factories_registration.hh
  LISTNAME   FLOW_ICY_OVERLAND_EVAL_FAC
  INSTALL    True
  )

