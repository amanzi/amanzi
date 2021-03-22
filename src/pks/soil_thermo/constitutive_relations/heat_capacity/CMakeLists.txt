# -*- mode: cmake -*-

#
#  ATS
#    Soil Thermo PK class
#

include_directories(${ATS_SOURCE_DIR}/src/factory)

add_library(soil_thermo_relations_heat_capacity
            soil_heat_capacity_evaluator.cc
            )

install(TARGETS soil_thermo_relations_heat_capacity DESTINATION lib)
