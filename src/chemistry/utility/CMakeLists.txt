# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Utilities
#

# Chemistry includes
include_directories(${CHEM_INCLUDES_DIR})

set(SRCS 
	 ./block.cc 
	 ./file_io.cc 
	 ./glenn.cc 
	 ./chemistry_exception.cc
	 ./string_tokenizer.cc
	 ./verbosity.cc
)

add_library(geochemutil ${SRCS}) 
add_install_library(geochemutil)


