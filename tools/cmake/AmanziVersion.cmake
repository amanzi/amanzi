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

message(STATUS ">>>>>>>> AmanziVersion.cmake")

find_package(Mercurial)
if ( (EXISTS ${CMAKE_SOURCE_DIR}/.hg/) AND (MERCURIAL_FOUND) ) 

  mercurial_branch(AMANZI_HG_BRANCH)
  mercurial_global_id(AMANZI_HG_GLOBAL_HASH)
  mercurial_local_id(AMANZI_HG_LOCAL_ID)

# set(MERCURIAL_ARGS head ${AMANZI_HG_BRANCH} --template "{latesttag('re:^amanzi-.*')}\n")
  set(MERCURIAL_ARGS log --rev "branch\(default\) and tag()" --template "{tags}\n" )
  execute_process(COMMAND  ${MERCURIAL_EXECUTABLE} ${MERCURIAL_ARGS}
	          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  RESULT_VARIABLE err_occurred 
                  OUTPUT_VARIABLE AMANZI_HG_LATEST_TAG
                  ERROR_VARIABLE err
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)

   # MESSAGE(STATUS ">>>> JDM: MERCURIAL_EXEC       = ${MERCURIAL_EXECUTABLE}")
   # MESSAGE(STATUS ">>>> JDM: MERCURIAL_ARGS       = ${MERCURIAL_ARGS}")
   # MESSAGE(STATUS ">>>> JDM: RESULT_VARIABLE      = ${err_occurred}")

   # Put the tags in a list
   STRING(REPLACE "\n" ";" AMANZI_HG_LATEST_TAG_LIST ${AMANZI_HG_LATEST_TAG})
   # Extract the lastest tag of the form amanzi-*
   FOREACH(atag ${AMANZI_HG_LATEST_TAG_LIST})
     IF ( ${atag} MATCHES "^amanzi-.*" )
       set ( AMANZI_HG_LATEST_TAG ${atag} )
     ENDIF()
   ENDFOREACH()

   STRING(REGEX REPLACE "amanzi-" "" AMANZI_HG_LATEST_TAG_VER ${AMANZI_HG_LATEST_TAG})	
   STRING(REGEX REPLACE "\\..*" "" AMANZI_HG_LATEST_TAG_MAJOR ${AMANZI_HG_LATEST_TAG_VER})	
   STRING(REGEX MATCH "\\.[0-9][0-9][\\.,-]" AMANZI_HG_LATEST_TAG_MINOR ${AMANZI_HG_LATEST_TAG_VER})	
   STRING(REGEX REPLACE "[\\.,-]" "" AMANZI_HG_LATEST_TAG_MINOR ${AMANZI_HG_LATEST_TAG_MINOR} )	

   set(AMANZI_VERSION_MAJOR ${AMANZI_HG_LATEST_TAG_MAJOR})
   set(AMANZI_VERSION_MINOR ${AMANZI_HG_LATEST_TAG_MINOR})

   #
   # Amanzi version
   #
   set(AMANZI_VERSION ${AMANZI_HG_LATEST_TAG_VER}_${AMANZI_HG_GLOBAL_HASH})

   STRING(REGEX REPLACE ".*\\.[0-9][0-9][\\.,-]" "" AMANZI_VERSION_PATCH ${AMANZI_VERSION})
   STRING(REGEX REPLACE ".*_" "" AMANZI_VERSION_HASH ${AMANZI_VERSION_PATCH})
   STRING(REGEX REPLACE "_.*" "" AMANZI_VERSION_PATCH ${AMANZI_VERSION_PATCH})

else()

  message(STATUS "  >>>>>>>> Using static version information to create amanzi_version.hh")
  if ( NOT MERCURIAL_FOUND ) 
    message(STATUS "    >>>>>> Could not locate Mercurial executable.")
  endif()
  if ( NOT EXISTS ${CMAKE_SOURCE_DIR}/.hg/ )
    message(STATUS "    >>>>>> Release or snapshot, no .hg directory found.")
  endif()

  #
  # For releases and snapshots, set static information before creating the tarball.
  #
  set(AMANZI_HG_BRANCH default )
  set(AMANZI_HG_GLOBAL_HASH )
  set(AMANZI_HG_LOCAL_ID )

  set(AMANZI_VERSION_MAJOR 0)
  set(AMANZI_VERSION_MINOR 85)
  set(AMANZI_VERSION_PATCH dev)
  set(AMANZI_VERSION_HASH ${AMANZI_HG_GLOBAL_HASH})

  #
  # Amanzi version
  #
  set(AMANZI_VERSION ${AMANZI_VERSION_MAJOR}.${AMANZI_VERSION_MINOR}-${AMANZI_VERSION_PATCH}_${AMANZI_VERSION_HASH})

endif()

# Write the version header filed
set(version_template ${AMANZI_SOURCE_TOOLS_DIR}/cmake/amanzi_version.hh.in)
configure_file(${version_template}
               ${CMAKE_CURRENT_BINARY_DIR}/amanzi_version.hh
               @ONLY)
configure_file(${version_template}
               ${CMAKE_CURRENT_BINARY_DIR}/extras/amanzi_version.hh
               @ONLY)

add_install_include_file(${CMAKE_CURRENT_BINARY_DIR}/amanzi_version.hh)             

message(STATUS "\t >>>>>  Amanzi Version: ${AMANZI_VERSION}")
message(STATUS "\t >>>>>  MAJOR ${AMANZI_VERSION_MAJOR}")
message(STATUS "\t >>>>>  MINOR ${AMANZI_VERSION_MINOR}")
message(STATUS "\t >>>>>  PATCH ${AMANZI_VERSION_PATCH}")
message(STATUS "\t >>>>>  HASH  ${AMANZI_VERSION_HASH}")

