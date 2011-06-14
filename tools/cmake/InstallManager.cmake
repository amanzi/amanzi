# # -*- mode: cmake -*-

#
# Functions for managing the install targets
#


include(CMakeParseArguments)
include(PrintVariable)

#
# Usage: ADD_INSTALL_INCLUDE_FILE( file1 file2 file3 ... )
#
# Arguments:
#  A list of files that will be installed in the AMANZI_INSTALL_INCLUDE_DIR
#
# 
function ( ADD_INSTALL_INCLUDE_FILE )

    foreach(_inc_file ${ARGV})
	install(FILES ${_inc_file} DESTINATION include)
    endforeach()	

endfunction( ADD_INSTALL_INCLUDE_FILE )

#
# Usage: ADD_INSTALL_LIBRARY( lib1 lib2 lib3 ... )
#
# Arguments:
#  A list of libraries that will be installed in the AMANZI_INSTALL_LIB_DIR
#
# 
function ( ADD_INSTALL_LIBRARY )

	install(TARGETS ${ARGV}
	        LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
		)

endfunction( ADD_INSTALL_LIBRARY )

#
# Usage: ADD_INSTALL_BINARY( exe1 exe2 ... )
#
# Arguments:
#  A list of executables that will be installed in the AMANZI_INSTALL_BIN_DIR
#
# 
function ( ADD_INSTALL_BINARY )
    
    foreach(_bin_file ${ARGV})
	#install(PROGRAMS ${_bin_file} DESTINATION bin)
	message(STATUS "BROKEN will nto install ${_bin_file}")
    endforeach()	

endfunction( ADD_INSTALL_BINARY )


