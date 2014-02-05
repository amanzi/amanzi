# -*- mode: cmake -*-

# ATS Surface balance PKs describe Evaporation, energy fluxes from
#  long/showtwave radiation, precip, etc etc etc


add_subdirectory(SEB)

#================================================
# register evaluators/factories/pks

# SEB

register_evaluator_with_factory(
  HEADERFILE SEB/surface_balance_SEB_reg.hh
  LISTNAME   SURFACE_BALANCE_SEB_REG
  )

generate_evaluators_registration_header(
  HEADERFILE surface_balance_SEB_registration.hh
  LISTNAME   SURFACE_BALANCE_SEB_REG
  INSTALL    True
  )

