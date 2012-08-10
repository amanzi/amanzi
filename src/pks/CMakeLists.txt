# -*- mode: cmake -*-

#
#  ATS
#    PK class
#

# ATS include directories
add_library(pk_factory pk_factory.cc)

add_subdirectory(mpc)
#add_subdirectory(energy)
#add_subdirectory(transport)
add_subdirectory(flow)




