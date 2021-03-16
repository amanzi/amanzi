# -*- mode: cmake -*-

#
#  ATS
#    Soil Thermo PK class
#

include_directories(${ATS_SOURCE_DIR}/src/factory)

add_library(soil_thermo_relations_thermal_conductivity
            soil_thermal_conductivity_evaluator.cc
            )

install(TARGETS soil_thermo_relations_thermal_conductivity DESTINATION lib)
