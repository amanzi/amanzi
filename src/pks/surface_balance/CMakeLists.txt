# -*- mode: cmake -*-

# ATS Surface balance PKs describe Evaporation, energy fluxes from
#  long/showtwave radiation, precip, etc etc etc

include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/surface_subsurface_fluxes)
include_directories(${ATS_SOURCE_DIR}/src/pks)

add_library(pk_surface_balance_SEB
       SEB/surface_balance_SEB.cc
       SEB/SnowEnergyBalance.cc
       SEB/surface_balance_evaluator.cc
       SEB_VPL/surface_balance_SEB_VPL.cc
       SEB_VPL/SnowEnergyBalance_VPL.cc
       SEB_VPL/surface_balance_evaluator_VPL.cc)

install(TARGETS pk_surface_balance_SEB DESTINATION lib)

#================================================
# register evaluators/factories/pks

# SEB
register_evaluator_with_factory(
  HEADERFILE SEB/surface_balance_SEB_reg.hh
  LISTNAME   SURFACE_BALANCE_SEB_REG
  )

# SEB VPL
register_evaluator_with_factory(
  HEADERFILE SEB_VPL/surface_balance_SEB_VPL_reg.hh
  LISTNAME   SURFACE_BALANCE_SEB_REG
  )

generate_evaluators_registration_header(
  HEADERFILE surface_balance_SEB_registration.hh
  LISTNAME   SURFACE_BALANCE_SEB_REG
  INSTALL    True
  )
