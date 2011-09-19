# # -*- mode: cmake -*-

#
# Functions for managing the install targets
#


include(CMakeParseArguments)
include(AmanziLinkLine)

export(PACKAGE Amanzi)

#
# Usage: ADD_INSTALL_INCLUDE_FILE( file1 file2 file3 ... )
#
# Arguments:
#  A list of files that will be installed in the AMANZI_INSTALL_INCLUDE_DIR
#
#
function ( ADD_INSTALL_INCLUDE_FILE )

  foreach(_inc_file ${ARGV})
    install(
      FILES ${_inc_file}
      DESTINATION include
      )
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

  install(
    TARGETS ${ARGV}
    EXPORT AmanziTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

  # Add the libraries to our global list
  add_amanzi_libraries(${ARGV})

  # Add dependency libaries as determined by the pacakge definition.
  add_package_libraries()

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
  install(
    TARGETS ${_bin_file}
    EXPORT AmanziTargets
    DESTINATION bin
    )
endforeach()

endfunction( ADD_INSTALL_BINARY )

#
# Usage: create_tpl_export_file( <package list> | package1 package2 ... )
#
# Arguments: Semicolon deliminated list of TPL package names or
#            individual package names
#
function( CREATE_TPL_EXPORT_FILE )

  # Print the usage
  macro( _print_usage )
    message("\nUsage: create_tpl_export_file( outfile PACKAGES <package list> | package1 package2 ... )\n")
  endmacro()

  # Parse Arguments
  set(_options "")
  set(_oneValue "")
  set(_multiValue "PACKAGES")
  cmake_parse_arguments(BUILD_TPL "${_options}" "${_oneValue}" "${_multiValue}" ${ARGN})

  if (NOT BUILD_TPL_PACKAGES)
    _print_usage()
    message(FATAL_ERROR "Require a package list to build export file")
  endif()

  list(GET BUILD_TPL_UNPARSED_ARGUMENTS 0 BUILD_TPL_OUTFILE)
  if (NOT BUILD_TPL_OUTFILE)
    _print_usage()
    message(FATAL_ERROR "Must define an output file")
  endif()  

  # BEGIN MACROS

  # Write the header for the file
  macro(_write_header)

    file(WRITE ${BUILD_TPL_OUTFILE} "# ############################################################################ #\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "#\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "# Amanzi TPL (External Software Packages) Configuration\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "#\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "# ############################################################################ #\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "\n")

  endmacro(_write_header)

  # Write a CMake set variable command
  macro(_write_cmake_variable _cmake_var_name _var_value)

    if (${_var_value})
      file(APPEND ${BUILD_TPL_OUTFILE} "set(${_cmake_var_name} ${${_var_value}})\n")
    endif()  

  endmacro(_write_cmake_variable)  

  # Add package to the file
  macro( _add_package _package )
	
    file(APPEND ${BUILD_TPL_OUTFILE} "\n#\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "# TPL: ${_package}\n")
    file(APPEND ${BUILD_TPL_OUTFILE} "#\n")

    set(found_package_flag ${_package}_FOUND)
    set(_package_enabled_flag_var "Amanzi_TPL_${_package}_ENABLED")
    set(_package_dir_var          "Amanzi_TPL_${_package}_DIR")

    if ( ${found_package_flag} ) 

      file(APPEND ${BUILD_TPL_OUTFILE} "set(${_package_enabled_flag_var} ON)\n")
      file(APPEND ${BUILD_TPL_OUTFILE} "set(${_package_dir_var} ${${_package}_DIR})\n")

      set(_var_name_list "INCLUDE_DIR;INCLUDE_DIRS;LIBRARIES;LIBRARY_DIR;LIBRARY_DIRS")
      foreach (_append_name ${_var_name_list})
	set(_write_var_name "Amanzi_TPL_${_package}_${_append_name}")
	set(_var_name       "${_package}_${_append_name}")
	_write_cmake_variable(${_write_var_name} ${_var_name})
      endforeach(_append_name)	

    else()

      file(APPEND ${BUILD_TPL_OUTFILE} "set(${_package_enabled_flag_var} OFF)\n")

    endif()

  endmacro(_add_package) 


  # Begin MAIN

  
  if (NOT EXISTS ${BUILD_TPL_OUTFILE})
    _write_header(${BUILD_TPL_OUTFILE})
  endif()  

  # Loop through each package and update the output file
  foreach(_package IN LISTS BUILD_TPL_PACKAGES)
    #print_variable(_package)
    _add_package(${_package}) 
  endforeach() 

  message(STATUS "Leaving create_tpl_export_file") 

endfunction ( CREATE_TPL_EXPORT_FILE )


#
# Usage: create_exports
#
# Arguments: None
#
#
function (CREATE_EXPORTS)

# Template file located in the CMake module directory

# Find the packages found for Amanzi
get_property(AMANZI_TPL_LIST GLOBAL PROPERTY PACKAGES_FOUND)
get_property(LINK_LINE GLOBAL PROPERTY AMANZI_LINK_LINE)

# Convert the link line to a space deliminated string
foreach (arg ${LINK_LINE})
  set(LINK_LINE_STRING "${LINK_LINE_STRING} ${arg}")
endforeach()

# Write and install the link-line file
file(WRITE ${AMANZI_LINK_LINE_FILE} ${LINK_LINE_STRING})
install(FILES ${AMANZI_LINK_LINE_FILE} DESTINATION lib)

# Write the TPL file
set(tpl_config_file "${AMANZI_BINARY_DIR}/AmanziConfigTPL.cmake")
if ( EXISTS ${tpl_config_file} )
  file(REMOVE ${tpl_config_file})
endif()  
create_tpl_export_file(${tpl_config_file}
                       PACKAGES ${AMANZI_ENABLED_TPLS})
install(FILES ${tpl_config_file} DESTINATION lib)				   

# Write the export Makefile and add to the include install list
set(in_makefile  "${AMANZI_MODULE_PATH}/MakefileConfig.export.in")
set(out_makefile "${AMANZI_BINARY_DIR}/Makefile.export")
configure_file("${in_makefile}" "${out_makefile}")
install(FILES "${out_makefile}" DESTINATION lib)

# Write the AmanziConfig.cmake file
set(in_config   "${AMANZI_MODULE_PATH}/AmanziConfig-install.cmake.in")
set(out_config   "${AMANZI_BINARY_DIR}/AmanziConfig.cmake")
configure_file(${in_config} ${out_config})
install(FILES ${out_config} DESTINATION lib)

# Write the AmanziConfigVersion.cmake file
set(in_config   "${AMANZI_MODULE_PATH}/AmanziConfigVersion-install.cmake.in")
set(out_config   "${AMANZI_BINARY_DIR}/AmanziConfigVersion.cmake")
configure_file(${in_config} ${out_config} @ONLY)
install(FILES ${out_config} DESTINATION lib)

# Write the CMake configuration target file
message(STATUS "Writing target file")
install(EXPORT AmanziTargets
        DESTINATION lib
	NAMESPACE amanzi_
	FILE AmanziTargets.cmake)

endfunction()
