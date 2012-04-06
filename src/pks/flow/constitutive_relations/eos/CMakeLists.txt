# -*- mode: cmake -*-

#
#  ATS
#    EOS for flow
#

include_directories(${ATS_SOURCE_DIR}/src/factory)
add_library(flow_relations_eos eos_factory.cc vapor_pressure_model_factory.cc eos_water.cc eos_ideal_gas.cc water_vapor_pressure_model.cc eos_vapor_in_gas.cc )


