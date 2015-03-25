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
#  snow_distribution_evaluator.cc
#  explicit_snow_distribution_evaluator.cc
#  implicit_snow_distribution_evaluator.cc
#  implicit_snow_distribution_uwflux_evaluator.cc
  snow_skin_potential_evaluator.cc
  )

install(TARGETS flow_relations_elevation DESTINATION lib)




