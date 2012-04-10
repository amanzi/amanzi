# -*- mode: cmake -*-

#
#  ATS
#    EOS for flow
#

include_directories(${ATS_SOURCE_DIR}/src/factory)
add_library(flow_relations_eos_ eos_factory.cc vapor_pressure_model_factory.cc eos_ideal_gas.cc eos_water.cc water_vapor_pressure_model.cc eos_vapor_in_gas.cc )
#add_library(flow_relations_eos_ eos_factory.cc eos_water.cc eos_ideal_gas.cc )

# this is very strange... our factory implementation relies on static
# variables which are initialized in cc files which then get included
# in archives.  These archives do not automatically include the needed
# symbols to initialize the BSS data, and therefore do not get sucked
# into the executable?  To get them, we need the whole archive...
#
# this is g++ specific and needs to be generalized!
set(flow_relations_eos -Wl,--whole-archive flow_relations_eos_ -Wl,--no-whole-archive)


if ( BUILD_TESTS )
    include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/eos)

    # Test: miscaleneous transport routines
    add_executable(eos_test test/main.cc)
    target_link_libraries(eos_test ${flow_relations_eos} ${Amanzi_TPL_Trilinos_LIBRARIES})
    #target_link_libraries(eos_test flow_relations_eos_factory ${Amanzi_TPL_Trilinos_LIBRARIES})

endif()
