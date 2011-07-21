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
	install(TARGETS ${_bin_file} DESTINATION bin)
    endforeach()	

endfunction( ADD_INSTALL_BINARY )

#
# Usage: create_export_makefile()
#
# Arguments: None
#
#
function (CREATE_EXPORT_MAKEFILE)

    # Template file locate in the CMake module directory
    set(in_makefile "${AMANZI_MODULE_PATH}/MakefileConfig.export.in")
    set(out_makefile "${AMANZI_BINARY_DIR}/Makefile.export")

    # Find the packages found for Amanzi
    get_property(AMANZI_TPL_LIST GLOBAL PROPERTY PACKAGES_FOUND) 

    # Write the export Makefile and add to the include install list
    configure_file("${in_makefile}" "${out_makefile}")
    install(FILES "${out_makefile}" DESTINATION lib)


endfunction(CREATE_EXPORT_MAKEFILE)
