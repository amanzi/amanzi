# -*- mode: cmake -*-

#
#  ATS
#    State class
#
message(STATUS "in state_dir: ATS_SOURCE_DIR=${ATS_SOURCE_DIR}")

# This state library will move to a new location
# once Markus has created the data manager.
# I use global properties since all the PROJECT_NAME
# variables in this directory are tied to the MPC.
# These properties will be replaced with the appropriate
# directory project name once it is moved.
#
# Library: state
#
add_library(state Field.cc State.cc Vis.cc)


