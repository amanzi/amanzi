# -*- mode: cmake -*-

#
#  ATS
#    WRM for flow
#


include_directories(${ATS_SOURCE_DIR}/src/factory)

add_library(flow_relations_wrm wrm_factory.cc wrm_linear_system.cc wrm_van_genuchten.cc)


