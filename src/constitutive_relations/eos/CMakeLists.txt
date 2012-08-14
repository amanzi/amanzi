# -*- mode: cmake -*-

#
#  ATS
#    EOS for flow
#

MACRO(WHOLE_ARCHIVE inlib)
  SET(${inlib} -Wl,--whole-archive ${inlib} -Wl,--no-whole-archive)
ENDMACRO(WHOLE_ARCHIVE inlib)

include_directories(${ATS_SOURCE_DIR}/src/factory)

add_library(flow_relations_eos
  eos_evaluator.cc
  eos_factory.cc
  eos_constant.cc eos_ideal_gas.cc eos_water.cc eos_ice.cc
  viscosity_evaluator.cc
  viscosity_relation_factory.cc
  viscosity_water.cc
  molar_fraction_gas_evaluator.cc
  vapor_pressure_relation_factory.cc
  vapor_pressure_water.cc
 )



#add_library(flow_relations_eos_factory eos_factory.cc)
#add_library(flow_relations_eos eos_water.cc eos_ideal_gas.cc )

# this is very strange... our factory implementation relies on static
# variables which are initialized in cc files which then get included
# in archives.  These archives do not automatically include the needed
# symbols to initialize the BSS data, and therefore do not get sucked
# into the executable?  To get them, we need the whole archive...
#
# this is g++ specific and needs to be generalized!
#WHOLE_ARCHIVE(flow_relations_eos_factory)
WHOLE_ARCHIVE(flow_relations_eos)

# if ( BUILD_TESTS )
#     include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/eos)

#     # Test: miscaleneous transport routines
#     add_executable(eos_test test/main.cc)
#     target_link_libraries(eos_test ${flow_relations_eos} ${Amanzi_TPL_Trilinos_LIBRARIES})
#     #target_link_libraries(eos_test ${flow_relations_eos_factory} ${flow_relations_eos} ${Amanzi_TPL_Trilinos_LIBRARIES})

# endif()
