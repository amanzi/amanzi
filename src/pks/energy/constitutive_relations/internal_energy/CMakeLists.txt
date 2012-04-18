# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#
#include_directories(${ATS_SOURCE_DIR}/src/factory)

add_library(energy_relations_internal_energy internal_energy_water_vapor.cc internal_energy_linear.cc)