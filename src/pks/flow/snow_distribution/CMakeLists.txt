# -*- mode: cmake -*-

include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/factory)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad/upwind_scheme)

add_library(
  pk_flow_snow_distribution
  snow_distribution_pk.cc 
  snow_distribution_ti.cc 
  snow_distribution_physics.cc
  )

install(TARGETS pk_flow_snow_distribution DESTINATION lib)

