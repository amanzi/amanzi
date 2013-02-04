# -*- mode: cmake -*-

#
#  ATS
#    WRM for flow
#

include_directories(${Amanzi_TPL_MSTK_INCLUDE_DIRS})
add_definitions("-DMSTK_HAVE_MPI")


add_library(flow_relations_elevation
  pres_elev_evaluator.cc
  elevation_evaluator.cc
  meshed_elevation_evaluator.cc
  standalone_elevation_evaluator.cc
  )




