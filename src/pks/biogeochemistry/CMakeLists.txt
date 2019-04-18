# -*- mode: cmake -*-

# ATS Surface balance PKs describe Evaporation, energy fluxes from
#  long/showtwave radiation, precip, etc etc etc


include_directories(${ATS_SOURCE_DIR}/src/factory)
include_directories(${ATS_SOURCE_DIR}/src/pks)

add_library(pk_BGC
  bgc_simple/utils.cc
  bgc_simple/PFT.cc
  bgc_simple/SoilCarbonParameters.cc
  bgc_simple/vegetation.cc
  bgc_simple/bgc_simple_funcs.cc
  bgc_simple/bgc_simple.cc
  carbon/simple/CarbonSimple.cc
  constitutive_models/carbon/bioturbation_evaluator.cc
)

install(TARGETS pk_BGC DESTINATION lib)

#================================================
# register evaluators/factories/pks

# SEB explicit
register_evaluator_with_factory(
  HEADERFILE bgc_simple/bgc_simple_reg.hh
  LISTNAME BGC_REG
)

register_evaluator_with_factory(
  HEADERFILE carbon/simple/CarbonSimple_reg.hh
  LISTNAME BGC_REG
)

register_evaluator_with_factory(
  HEADERFILE constitutive_models/carbon/bioturbation_evaluator_reg.hh
  LISTNAME BGC_REG
)

generate_evaluators_registration_header(
  HEADERFILE BGC_registration.hh
  LISTNAME   BGC_REG
  INSTALL    True
)


