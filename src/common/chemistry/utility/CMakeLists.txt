# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Utilities
#

# Chemistry includes
include_directories(${CHEM_INCLUDES_DIR})

set(SRCS 
	 ./matrix_block.cc 
	 ./chemistry_exception.cc
	 ./chemistry_utilities.cc
	 ./string_tokenizer.cc
	 ./chemistry_verbosity.cc
)

add_amanzi_library(geochemutil SOURCE ${SRCS}
	                       LINK_LIBS ${Teuchos_LIBRARIES}) 


