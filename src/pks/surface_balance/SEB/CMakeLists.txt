# -*- mode: cmake -*-

include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/surface_subsurface_fluxes)
include_directories(${ATS_SOURCE_DIR}/src/pks)
add_library(pk_surface_balance_SEB surface_balance_SEB.cc SnowEnergyBalance.cc)


