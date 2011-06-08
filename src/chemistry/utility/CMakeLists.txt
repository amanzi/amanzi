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
	 ./chemistry-exception.cc
	 ./string-tokenizer.cc
)

add_library(geochemutil ${SRCS}) 

