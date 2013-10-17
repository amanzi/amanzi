# -*- mode: cmake -*-

#
#  ATS
#

add_library(flow_relations_surface_subsurface_fluxes
  overland_source_from_subsurface_flux_evaluator.cc
  surface_top_cells_evaluator.cc
  )



register_evaluator_with_factory(
  HEADERFILE surface_top_cells_evaluator_reg.hh
  LISTNAME   CONSTITUTIVE_RELATIONS_SURFACE_SUBSURFACE_FLUXES_EVALUATORS
  )

register_evaluator_with_factory(
  HEADERFILE overland_source_from_subsurface_flux_evaluator_reg.hh
  LISTNAME   CONSTITUTIVE_RELATIONS_SURFACE_SUBSURFACE_FLUXES_EVALUATORS
  )

generate_evaluators_registration_header(
  HEADERFILE constitutive_relations_surface_subsurface_fluxes_evaluators_registration.hh
  LISTNAME   CONSTITUTIVE_RELATIONS_SURFACE_SUBSURFACE_FLUXES_EVALUATORS
  INSTALL    True
  )



