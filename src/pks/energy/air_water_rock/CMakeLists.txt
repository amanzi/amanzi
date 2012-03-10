# -*- mode: cmake -*-

include_directories(${ATS_SOURCE_DIR}/src/operators/advection)
add_library(energy_air_water_rock air_water_rock_pk.cc air_water_rock_ti.cc air_water_rock_physics.cc)

