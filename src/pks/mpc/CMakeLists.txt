# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

# ATS include directories
include_directories(${ATS_SOURCE_DIR}/src/data_structures)
include_directories(${ATS_SOURCE_DIR}/src/state)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/eos)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/ewc)
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/wrm)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/wrm/models)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/richards)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/overland)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/overland_head)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/constitutive_relations/internal_energy)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/base)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/two_phase)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad)

add_library(mpc
  weak_mpc.cc
  mpc_coupled_cells.cc
  mpc_surface_subsurface_coupler.cc
  mpc_surface_subsurface_dirichlet_coupler.cc
  mpc_surface_subsurface_flux_coupler.cc
  mpc_surface_subsurface_energy_coupler.cc
  mpc_delegate_ewc.cc
  mpc_subsurface.cc
  mpc_permafrost.cc)


# register factories
register_evaluator_with_factory(
  HEADERFILE weak_mpc_reg.hh
  LISTNAME   MPC_PK_FAC
  )

register_evaluator_with_factory(
  HEADERFILE strong_mpc_reg.hh
  LISTNAME   MPC_PK_FAC
  )

register_evaluator_with_factory(
  HEADERFILE mpc_subsurface_reg.hh
  LISTNAME   MPC_PK_FAC
  )

register_evaluator_with_factory(
  HEADERFILE mpc_permafrost_reg.hh
  LISTNAME   MPC_PK_FAC
  )

register_evaluator_with_factory(
  HEADERFILE mpc_surface_subsurface_dirichlet_coupler_reg.hh
  LISTNAME   MPC_PK_FAC
  )

register_evaluator_with_factory(
  HEADERFILE mpc_surface_subsurface_flux_coupler_reg.hh
  LISTNAME   MPC_PK_FAC
  )

register_evaluator_with_factory(
  HEADERFILE mpc_surface_subsurface_energy_coupler_reg.hh
  LISTNAME   MPC_PK_FAC
  )

register_evaluator_with_factory(
  HEADERFILE mpc_water_coupler_reg.hh
  LISTNAME   MPC_PK_FAC
  )

generate_evaluators_registration_header(
  HEADERFILE mpc_registration.hh
  LISTNAME   MPC_PK_FAC
  INSTALL    True
  )
