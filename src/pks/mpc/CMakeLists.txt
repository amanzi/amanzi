# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

# ATS include directories
include_directories(${ATS_SOURCE_DIR}/src/data_structures)
include_directories(${ATS_SOURCE_DIR}/src/state)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/eos)
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

add_library(mpc mpc.cc weak_mpc.cc strong_mpc.cc
  mpc_diagonal_flow_energy.cc
  mpc_coupled_flow_energy.cc
  mpc_frozen_coupled_flow_energy.cc
  permafrost_model.cc
  mpc_surface_subsurface_water_coupler.cc
  mpc_surface_subsurface_neumann_coupler.cc
  mpc_surface_subsurface_dirichlet_coupler.cc)


