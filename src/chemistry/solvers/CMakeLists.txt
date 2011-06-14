# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Reactions
#

# Chemistry includes
include_directories(${CHEM_INCLUDES_DIR})

set(SRCS 
	 ./LU.cpp 
	 ./Newton.cpp 
)

add_library(geochemsolvers ${SRCS}) 
add_install_library(geochemsolvers)

