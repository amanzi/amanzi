# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Reactions
#

# Chemistry includes
include_directories(${CHEM_INCLUDES_DIR})

set(geochemsolvers_src newton.cc)

set(geochemsolvers_inc lu_solver.hh
	               newton.hh)

add_amanzi_library(geochemsolvers SOURCE ${geochemsolvers_src}
	                          LINK_LIBS geochemutil) 

add_install_include_file(${geochemsolvers_inc})
