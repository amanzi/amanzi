# -*- mode: cmake -*-

#
#  ATS
#    WRM for flow
#

#MACRO(WHOLE_ARCHIVE inlib)
#  SET(${inlib} -Wl,--whole-archive ${inlib} -Wl,--no-whole-archive)
#ENDMACRO(WHOLE_ARCHIVE inlib)

#include_directories(${ATS_SOURCE_DIR}/src/factory)

add_library(flow_relations_pc pc_ice_water.cc)


