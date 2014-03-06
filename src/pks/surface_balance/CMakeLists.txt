# -*- mode: cmake -*-

# ATS Surface balance PKs describe Evaporation, energy fluxes from
#  long/showtwave radiation, precip, etc etc etc


add_subdirectory(SEB)
add_subdirectory(SEB_VPL)


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


# SEB_VPL
register_evaluator_with_factory(
  HEADERFILE SEB_VPL/surface_balance_SEB_VPL_reg.hh
  LISTNAME   SURFACE_BALANCE_SEB_VPL_REG
  )

generate_evaluators_registration_header(
  HEADERFILE surface_balance_SEB_VPL_registration.hh
  LISTNAME   SURFACE_BALANCE_SEB_VPL_REG
  INSTALL    True
  )



