# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Utilities
#

# Chemistry includes
include_directories(${CHEM_INCLUDES_DIR})

set(SRCS 
	 ./Block.cpp 
	 ./FileIO.cpp 
	 ./Glenn.cpp 
	 ./ChemistryException.cpp
	 ./StringTokenizer.cpp
)

add_library(geochemutil ${SRCS}) 

