# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

include_directories(${ATS_SOURCE_DIR}/src/factory)

add_library(energy_relations_thermal_conductivity thermal_conductivity_factory.cc thermal_conductivity_twophase_peterslidard.cc thermal_conductivity_twophase_wetdry.cc)