# -*- mode: cmake -*-

#
#  ATS
#    internal energy constitutive relations
#
include_directories(${ATS_SOURCE_DIR}/src/factory)

add_library(energy_relations_internal_energy
             iem_evaluator.cc
             iem_factory.cc
             iem_linear.cc
             iem_quadratic.cc
             iem_water_vapor_evaluator.cc
             iem_water_vapor.cc
             )

install(TARGETS energy_relations_internal_energy DESTINATION lib)

