# -*- mode: cmake -*-

include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/surface_subsurface_fluxes)
include_directories(${ATS_SOURCE_DIR}/src/pks)

add_library(pk_surface_balance_SEB_VPL surface_balance_SEB_VPL.cc SnowEnergyBalance_VPL.cc surface_balance_evaluator_VPL.cc)

install(TARGETS pk_surface_balance_SEB_VPL DESTINATION lib)


