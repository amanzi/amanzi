# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Reactions
#

# Chemistry includes
include_directories(${CHEM_INCLUDES_DIR})

set(SRCS 
	 ./lu_solver.cc
	 ./newton.cc 
)

add_amanzi_library(geochemsolvers SOURCE ${SRCS}) 
