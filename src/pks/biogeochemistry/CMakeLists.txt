# -*- mode: cmake -*-

# ATS Surface balance PKs describe Evaporation, energy fluxes from
#  long/showtwave radiation, precip, etc etc etc


include_directories(${ATS_SOURCE_DIR}/src/pks)

add_library(pk_biogeochemistry
  bgc_simple/utils.cc
  bgc_simple/SoilCarbonParameters.cc
  bgc_simple/vegetation.cc
  bgc_simple/bgc_simple_funcs.cc
)

install(TARGETS pk_biogeochemistry DESTINATION lib)

