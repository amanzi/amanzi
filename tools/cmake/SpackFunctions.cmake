# -*- mode: cmake -*-

#
# Amanzi Spack Information:
# 
# Used for leveraging Spack functionality to install packages, find 
# install locations for Spack packages and their version information
#
# 
#


message(STATUS ">>>>>>>> SpackFunctions.cmake")

#########################################################
#
# Install the package
#
#########################################################

function(spack_install_package pkg)
  message(STATUS "   >>>>>>>> Installing ${pkg}")
  
  set(SPACK_ARGS install ${pkg})
  execute_process(COMMAND  ${SPACK_BINARY} ${SPACK_ARGS}
	          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  RESULT_VARIABLE err_occurred 
                  OUTPUT_VARIABLE spack_status
                  ERROR_VARIABLE err
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)
  unset(SPACK_ARGS)

  if(err_occurred)
    message(WARNING "Error executing spack install:\n ${cmd}\n${err}")
    exit()
  endif()
endfunction(spack_install_package)

#########################################################
#
# Create symlinks between where spack installed the 
# package and where you actually want it installed
# (i.e. the directory specified by TPL_INSTALL_PREFIX)
#
#########################################################

function ( spack_create_symlinks pkg install_prefix )
  message(STATUS "   >>>>>>>> Making symlinks for ${pkg}")
  
  set(SPACK_ARGS view symlink ${install_prefix} ${pkg})
  execute_process(COMMAND  ${SPACK_BINARY} ${SPACK_ARGS}
	          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  RESULT_VARIABLE err_occurred 
                  OUTPUT_VARIABLE spack_status
                  ERROR_VARIABLE err
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)
  unset(SPACK_ARGS)

  if(err_occurred)
    message(WARNING "Error executing spack view:\n ${cmd}\n${err}")
    exit()
  endif()
endfunction(spack_create_symlinks)

#########################################################
#
# Alternatively, locate where Spack installed the package
# instead of creating symlinks within the desired install
# location...
#
# This sets the variable TMP_SPACK_INSTALL_DIR to the full
# path of the install directory of the given package. 
# 
# Afterwards, variables of the form <package>_INSTALL_PREFIX
# and <package>_INCLUDE_DIRS (and so on) can be set using 
# the info contained in this variable.
#
########################################################
function ( spack_locate_pkg_install_dir pkg )
  message(STATUS "   >>>>>>>> Locating where spack installed ${pkg}")
  
  set(SPACK_ARGS location -i ${pkg})
  execute_process(COMMAND  ${SPACK} ${SPACK_ARGS}
	          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  RESULT_VARIABLE err_occurred 
                  OUTPUT_VARIABLE TMP_SPACK_INSTALL_DIR
                  ERROR_VARIABLE err
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)
  unset(SPACK_ARGS)

  if(err_occurred)
    message(WARNING "Error executing spack location:\n ${cmd}\n${err}")
    exit()
  endif()

  set(TMP_SPACK_INSTALL_DIR ${TMP_SPACK_INSTALL_DIR} PARENT_SCOPE)

  message(STATUS "  >>>>>>>> Directory found: ${TMP_SPACK_INSTALL_DIR}")
endfunction(spack_locate_pkg_install_dir)

########################################################
#
# Extract what version of the package was installed
#
# This sets the variable TMP_SPACK_VERSION to the version
# of the package that was installed
# 
# Afterwards, variables of the form <package>_VERSION_MAJOR
# and <package>_VERSION_MINOR (and so on) can be set using
# the info contained in this variable.
#
#########################################################
function ( spack_extract_version_info pkg )
  message(STATUS "   >>>>>>>> extracting version information for ${pkg}")
  
  set(SPACK_ARGS find ${pkg})
  execute_process(COMMAND  ${SPACK} ${SPACK_ARGS}
	          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  RESULT_VARIABLE err_occurred 
                  OUTPUT_VARIABLE TMP_SPACK_VERSION
                  ERROR_VARIABLE err
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)
  unset(SPACK_ARGS)

  if(err_occurred)
    message(WARNING "Error executing spack find:\n ${cmd}\n${err}")
    exit()
  endif()

  STRING(REGEX REPLACE "^${pkg}@" "" TMP_SPACK_VERSION ${TMP_SPACK_VERSION})

  set(TMP_SPACK_VERSION ${TMP_SPACK_VERSION} PARENT_SCOPE )
  message(STATUS "  >>>>>>>> Version found: ${TMP_SPACK_VERSION}")
endfunction(spack_extract_version_info)
