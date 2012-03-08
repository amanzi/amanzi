# -*- mode: cmake -*-

#
#  ATS
#    MPC class
#
# ATS include directories
include_directories(${ATS_SOURCE_DIR}/src/data_structures)
include_directories(${ATS_SOURCE_DIR}/src/state)
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow)


add_library(coordinator Coordinator.cc Unstructured_observations.cc)

