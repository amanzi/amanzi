# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

# ATS include directories
include_directories(${ATS_SOURCE_DIR}/src/state)
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/utils/eos)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/utils/saturation_curves)

# This state library will move to a new location
# once Markus has created the data manager.
# I use global properties since all the PROJECT_NAME
# variables in this directory are tied to the MPC.
# These properties will be replaced with the appropriate
# directory project name once it is moved.
#
# Library: state
#
add_library(permafrost PermafrostProblem.cc PermafrostModelEvaluator.cc Permafrost_PK.cc)


