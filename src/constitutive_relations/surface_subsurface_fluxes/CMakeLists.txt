# -*- mode: cmake -*-

#
#  ATS
#

add_library(flow_relations_surface_subsurface_fluxes
  overland_source_from_subsurface_flux_evaluator.cc
  surface_top_cells_evaluator.cc
  top_cells_surface_evaluator.cc
  volumetric_darcy_flux_evaluator.cc
  )

install(TARGETS flow_relations_surface_subsurface_fluxes DESTINATION lib )


