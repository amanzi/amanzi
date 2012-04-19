# -*- mode: cmake -*-

#
#  ATS
#    WRM for flow
#

MACRO(WHOLE_ARCHIVE inlib)
  SET(${inlib} -Wl,--whole-archive ${inlib} -Wl,--no-whole-archive)
ENDMACRO(WHOLE_ARCHIVE inlib)

include_directories(${ATS_SOURCE_DIR}/src/factory)

add_library(flow_relations_wrm wrm_factory.cc wrm_linear_system.cc wrm_van_genuchten.cc)
WHOLE_ARCHIVE(flow_relations_wrm)

if ( BUILD_TESTS )
    include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/wrm)

    # Test: miscaleneous transport routines
    add_executable(plot_wrm test/plot_wrm.cc)
    target_link_libraries(plot_wrm ${flow_relations_wrm} ${Amanzi_TPL_Trilinos_LIBRARIES})
endif()

