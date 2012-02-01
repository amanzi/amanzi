# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

# ATS include directories
include_directories(${ATS_SOURCE_DIR}/src/state)
include_directories(${ATS_SOURCE_DIR}/src/data_structures)
include_directories(${ATS_SOURCE_DIR}/src/time_integration)
include_directories(${ATS_SOURCE_DIR}/src/ats_solvers)
include_directories(${ATS_SOURCE_DIR}/src/pks)

add_library(energy NullEnergyPK.cc)


