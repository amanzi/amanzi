# -*- mode: cmake -*-

#
#  ATS
#    internal energy constitutive relations
#
include_directories(${ATS_SOURCE_DIR}/src/factory)

MACRO(WHOLE_ARCHIVE inlib)
  SET(${inlib} -Wl,--whole-archive ${inlib} -Wl,--no-whole-archive)
ENDMACRO(WHOLE_ARCHIVE inlib)

add_library(energy_relations_internal_energy 
             internal_energy_water_vapor.cc 
             iem_factory.cc 
             iem_linear.cc 
             iem_quadratic.cc
             )