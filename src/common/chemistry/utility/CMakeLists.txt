# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Utilities
#

# Chemistry includes
include_directories(${CHEM_INCLUDES_DIR})

set(geochemutil_src_files matrix_block.cc 
	                  chemistry_utilities.cc
)

set(geochemutil_inc_files matrix_block.hh
                          chemistry_exception.hh
                          chemistry_utilities.hh)

add_amanzi_library(geochemutil SOURCE ${geochemutil_src_files}
	                       LINK_LIBS ${Teuchos_LIBRARIES}) 

# set_target_properties(geochemutil PROPERTIES PUBLIC_HEADER ${geochemutil_inc_files})

add_install_include_file(${geochemutil_inc_files})

