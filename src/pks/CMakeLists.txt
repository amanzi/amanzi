# -*- mode: cmake -*-

#
#  ATS
#    PK class
#

#add_subdirectory(mpc)
add_subdirectory(energy)
#add_subdirectory(flow)

# ATS include directories
include_directories(${ATS_SOURCE_DIR}/src/pks/mpc)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/constant_temperature)
include_directories(${ATS_SOURCE_DIR}/src/operators/advection)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/air_water_rock)

add_library(pks PK_Factory.cc)



