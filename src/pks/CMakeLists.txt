# -*- mode: cmake -*-

#
#  ATS
#    PK class
#

# ATS include directories
add_library(pk_bases
  pk_physical_bdf_base.cc
  pk_physical_base.cc
  pk_bdf_base.cc
  pk_default_base.cc
  pk_factory.cc)

add_subdirectory(mpc)
add_subdirectory(flow)
add_subdirectory(energy)
#add_subdirectory(transport)
add_subdirectory(deform)





