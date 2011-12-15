# -*- mode: cmake -*-

#
#  ATS
#    PK class
#

add_subdirectory(flow)

# ATS include directories
include_directories(${ATS_SOURCE_DIR}/src/state)
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/Darcy)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/Richards)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/permafrost)

# This state library will move to a new location
# once Markus has created the data manager.
# I use global properties since all the PROJECT_NAME
# variables in this directory are tied to the MPC.
# These properties will be replaced with the appropriate
# directory project name once it is moved.
#
# Library: state
#
add_library(pks PK_Factory.cc)



