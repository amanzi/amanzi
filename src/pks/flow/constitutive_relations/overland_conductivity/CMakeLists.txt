# -*- mode: cmake -*-

#
#  ATS
#    WRM for flow
#

MACRO(WHOLE_ARCHIVE inlib)
  SET(${inlib} -Wl,--whole-archive ${inlib} -Wl,--no-whole-archive)
ENDMACRO(WHOLE_ARCHIVE inlib)

add_library(flow_relations_overland_conductivity
  manning_conductivity_model.cc
  overland_conductivity_evaluator.cc
  unfrozen_fraction_model.cc
  unfrozen_fraction_evaluator.cc
  unfrozen_effective_depth_evaluator.cc
  )




