# -*- mode: cmake -*-

#
# Amanzi Version Information:
# 
# Information about the current source is extracted from the mercurial repository and used to 
# create the version string (AMANZI_VERSION).  
#
# NOTE: this information won't be accessible without the full repository.
#       So for releases we need to extract this and set it as part of the tarball creation.
#
#   * if amanzi_version.hh does not exist create it
#       * if mercurial is found
#            use mercurial to create version strings 
#       * else
#            use statically defined version strings
#       * endif
#   * endif
#   install amanzi_version.hh
#

include(MercurialMacros)
include(PrintVariable)
include(InstallManager)

find_package(Mercurial)
if ( NOT MERCURIAL_FOUND ) 

  message(ERROR "Could not locate Mercurial executable. Setting static version information")

  #
  # Not sure what is best for static information.  
  #
  set(AMANZI_VERSION_MAJOR 0)
  set(AMANZI_VERSION_MINOR 80-dev)
  set(AMANZI_VERSION_PATCH )

  #
  # Amanzi version
  #
  set(AMANZI_VERSION ${AMANZI_VERSION_MAJOR}.${AMANZI_VERSION_MINOR})

else()

  mercurial_branch(AMANZI_HG_BRANCH)
  mercurial_global_id(AMANZI_HG_GLOBAL_HASH)
  mercurial_local_id(AMANZI_HG_LOCAL_ID)

  set(MERCURIAL_ARGS head ${AMANZI_HG_BRANCH} --template {latesttag}\n)
  execute_process(COMMAND  ${MERCURIAL_EXECUTABLE} ${MERCURIAL_ARGS}
	          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  RESULT_VARIABLE err_occurred 
                  OUTPUT_VARIABLE AMANZI_HG_LATEST_TAG
                  ERROR_VARIABLE err
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)

   STRING(REGEX REPLACE "amanzi-" "" AMANZI_HG_LATEST_TAG_VER ${AMANZI_HG_LATEST_TAG})	

   #
   # Amanzi version
   #
   set(AMANZI_VERSION ${AMANZI_HG_LATEST_TAG_VER}_${AMANZI_HG_GLOBAL_HASH})

endif()

# Status output
message(STATUS "Amanzi Version:\t${AMANZI_VERSION}")
message(STATUS "\tMercurial Branch:\t${AMANZI_HG_BRANCH}")
message(STATUS "\tMercurial Global ID:\t${AMANZI_HG_GLOBAL_HASH}")
message(STATUS "\tMercurial Local ID:\t${AMANZI_HG_LOCAL_ID}")
message(STATUS "\tMercurial Tag Version:\t${AMANZI_HG_LATEST_TAG_VER}")


# Write the version header filed
set(version_template ${AMANZI_SOURCE_TOOLS_DIR}/cmake/amanzi_version.hh.in)
configure_file(${version_template}
               ${CMAKE_CURRENT_BINARY_DIR}/amanzi_version.hh
               @ONLY)

add_install_include_file(${CMAKE_CURRENT_BINARY_DIR}/amanzi_version.hh)             


