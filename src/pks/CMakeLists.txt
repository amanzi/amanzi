# -*- mode: cmake -*-

#
#  ATS
#    PK class
#

add_subdirectory(mpc)
add_subdirectory(energy)
#add_subdirectory(flow)

# ATS include directories
include_directories(${ATS_SOURCE_DIR}/src/data_structures)
include_directories(${ATS_SOURCE_DIR}/src/time_integration)
include_directories(${ATS_SOURCE_DIR}/src/ats_solvers)
include_directories(${ATS_SOURCE_DIR}/src/state)
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/mpc)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy)
#include_directories(${ATS_SOURCE_DIR}/src/pks/flow/Darcy)
#include_directories(${ATS_SOURCE_DIR}/src/pks/flow/Richards)
#include_directories(${ATS_SOURCE_DIR}/src/pks/flow/permafrost)

add_library(pks PK_Factory.cc)



