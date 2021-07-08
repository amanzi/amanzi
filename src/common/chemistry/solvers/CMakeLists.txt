# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Reactions
#

# Chemistry includes
include_directories(${CHEM_INCLUDES_DIR})

set(geochemsolvers_src Newton.cc)

set(geochemsolvers_inc LUSolver.hh
	               Newton.hh)

add_amanzi_library(geochemsolvers SOURCE ${geochemsolvers_src}
	                          LINK_LIBS geochemutil) 

add_install_include_file(${geochemsolvers_inc})
