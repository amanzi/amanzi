# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

# ATS include directories
add_subdirectory(passive_tracer)

#================================================
# register evaluators/factories/pks

# passive_tracer

register_evaluator_with_factory(
  HEADERFILE passive_tracer/passive_tracer_reg.hh
  LISTNAME   TRANSPORT_PASSIVE_TRACER_REG
  )

generate_evaluators_registration_header(
  HEADERFILE transport_passive_tracer_registration.hh
  LISTNAME   TRANSPORT_PASSIVE_TRACER_REG
  INSTALL    True
  )



