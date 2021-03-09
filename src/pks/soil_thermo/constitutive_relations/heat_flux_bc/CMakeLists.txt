# -*- mode: cmake -*-

#
#  ATS
#    Soil Thermo Heat Flux Evaluator class
#

include_directories(${ATS_SOURCE_DIR}/src/factory)

add_library(soil_thermo_relations_heat_flux_bc
            heat_flux_bc_evaluator.cc
            )

install(TARGETS soil_thermo_relations_heat_flux_bc DESTINATION lib)
