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

  # MACROS

  # Print the usage
  macro( _print_usage )
    message("\nUsage: create_tpl_export_file(<package list> | package1 package2 ... )\n")
  endmacro()

  # Write the header for the file
  macro(_write_header _outfile)

    file(WRITE ${_outfile} "# ############################################################################ #\n")
    file(APPEND ${_outfile} "#\n")
    file(APPEND ${_outfile} "# Amanzi TPL (External Software Packages) Configuration\n")
    file(APPEND ${_outfile} "#\n")
    file(APPEND ${_outfile} "# ############################################################################ #\n")
    file(APPEND ${_outfile} "\n")

  endmacro(_write_header)

  macro(_write_variable _var_name _var_value _outfile)

    file(APPEND ${_outfile} "set(${_var_name} ${_var_value})\n")

  endmacro(_write_variable)  

  # Add packge to the file
  macro( _add_package _package _outfile )

    set(found_package_flag ${_package}_FOUND)

    if ( ${found_package_flag} ) 
      set(_package_enabled_flag_var "Amanzi_TPL_${_package}_ENABLED")
      set(_package_dir_var          "Amanzi_TPL_${_package}_DIR")
      set(_package_include_dir_var  "Amanzi_TPL_${_package}_INCLUDE_DIR")
      set(_package_include_dirs_var "Amanzi_TPL_${_package}_INCLUDE_DIRS")
      set(_package_libraries_var    "Amanzi_TPL_${_package}_LIBRARIES")
      set(_package_library_dir_var  "Amanzi_TPL_${_package}_LIBRARY_DIR")
      set(_package_library_dirs_var "Amanzi_TPL_${_package}_LIBRARY_DIRS")

      file(APPEND ${_outfile} "\n#\n")
      file(APPEND ${_outfile} "# TPL: ${_package}\n")
      file(APPEND ${_outfile} "#\n")
      
      file(APPEND ${_outfile} "set(${_package_enabled_flag_var} ON)\n")
      file(APPEND ${_outfile} "set(${_package_dir_var}          ${${_package}_DIR})\n")

      # Include directories
      if ( ${_package}_INCLUDE_DIR ) 
	_write_variable(${_package_include_dir_var} ${${_package}_INCLUDE_DIR} ${_outfile})
      endif()	
      #file(APPEND ${_outfile} "set(${_package_include_dir_var}  ${${_package}_INCLUDE_DIR})\n")
      file(APPEND ${_outfile} "set(${_package_include_dirs_var} ${${_package}_INCLUDE_DIRS})\n")
     
      # Libraries and library directories
      file(APPEND ${_outfile} "set(${_package_libraries_var} ${${_package}_LIBRARIES})\n")
      if ( ${_package}_LIBRARY_DIR )
	_write_variable(${_package_library_dir_var} ${${_package}_LIBRARY_DIR} ${_outfile})
      endif()	
      if ( ${_package}_LIBRARY_DIRS )
	_write_variable(${_package_library_dirs_var} ${${_package}_LIBRARY_DIRS} ${_outfile})
      endif()	
      
    else()
      file(APPEND ${_outfile} "set(${_package_enabled_flag_var} OFF)\n")
    endif()

  endmacro(_add_package) 


  # Begin MAIN
  set(_package_list ${ARGV})

  if(NOT _package_list)
    message(FATAL_ERROR "\nUsage create_tpl_export_file(<package_list>| package1 package2..)\n")
  endif(NOT _package_list)

  set(outfile "${AMANZI_BINARY_DIR}/AmanziConfigTPL.cmake")
  install(FILES ${outfile} DESTINATION lib) 
  #print_variable(outfile)

  if (NOT EXISTS ${outfile})
    _write_header(${outfile})
  endif()  

  # Write the header
  _write_header(${outfile})

  foreach(_package IN LISTS _package_list)

    #print_variable(_package)
    _add_package(${_package} ${outfile}) 

  endforeach()  

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
create_tpl_export_file(${AMANZI_ENABLED_TPLS})

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
