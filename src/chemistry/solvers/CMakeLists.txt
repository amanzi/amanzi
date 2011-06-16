# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Reactions
#

# Chemistry includes
include_directories(${CHEM_INCLUDES_DIR})

set(SRCS 
	 ./lu.cc
	 ./newton.cc 
)

add_library(geochemsolvers ${SRCS}) 
add_install_library(geochemsolvers)

