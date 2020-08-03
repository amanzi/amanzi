# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

include_directories(${ATS_SOURCE_DIR}/src/factory)

add_library(energy_relations_thermal_conductivity
            thermal_conductivity_twophase_evaluator.cc
            thermal_conductivity_twophase_factory.cc
            thermal_conductivity_twophase_peterslidard.cc
            thermal_conductivity_twophase_wetdry.cc
            thermal_conductivity_threephase_evaluator.cc
            thermal_conductivity_threephase_factory.cc
            thermal_conductivity_threephase_peterslidard.cc
            thermal_conductivity_threephase_wetdry.cc
            thermal_conductivity_threephase_volume_averaged.cc
            thermal_conductivity_threephase_sutra_hacked.cc
            thermal_conductivity_surface_evaluator.cc
            )

install(TARGETS energy_relations_thermal_conductivity DESTINATION lib)
