# -*- mode: cmake -*-
include_directories(${ATS_SOURCE_DIR}/src/pks)
add_library(pk_energy_constant_temperature constant_temperature.cc)

