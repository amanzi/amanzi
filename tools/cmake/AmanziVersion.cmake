# -*- mode: cmake -*-

#
# Amanzi Version Information:
#
# Information about the current source is extracted from the git repository and used to
# create the version string (AMANZI_VERSION).
#
# NOTE: this information won't be accessible without the full repository.
#       So for releases we need to extract this and set it as part of the tarball creation.
#
#   * if amanzi_version.hh does not exist create it
#       * if git is found
#            use git to create version strings
#       * else
#            use statically defined version strings
#       * endif
#   * endif
#   install amanzi_version.hh
#

include(PrintVariable)
include(InstallManager)

message(STATUS "")
message(STATUS ">>>>>>>> AmanziVersion.cmake")

find_package(Git)

if ( (EXISTS ${CMAKE_SOURCE_DIR}/.git/) AND (GIT_FOUND) )

  # Get the name of the current branch.
  set(GIT_ARGS status)
  execute_process(COMMAND ${GIT_EXECUTABLE} ${GIT_ARGS}
	          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  RESULT_VARIABLE err_occurred
                  OUTPUT_VARIABLE AMANZI_GIT_STATUS
                  ERROR_VARIABLE err
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)
  if (err_occurred)
    message(WARNING "Error executing git:\n ${cmd}\n${err}")
    set(cmd_output cmd_output-NOTFOUND)
    exit()
  endif()

  # too noisy command...
  # message(STATUS ">>>> JDM: AMANZI_GIT_STATUS:      ${AMANZI_GIT_STATUS}")

  # Put the status in a list
  STRING(REPLACE "\n" ";" AMANZI_GIT_STATUS_LIST ${AMANZI_GIT_STATUS})
  # Extract the first entry - reuse the AMANZI_GIT_STATUS variable
  LIST(GET AMANZI_GIT_STATUS_LIST 0 AMANZI_GIT_STATUS)
  if (${AMANZI_GIT_STATUS} MATCHES "(D|d)etached")
    # For now just set branch to detached - we could add a lookup for tags later
    set(AMANZI_GIT_BRANCH detached)
  elseif(${AMANZI_GIT_STATUS} MATCHES "On branch")
    # Extract the branch name
    STRING(REPLACE "On branch " "" AMANZI_GIT_BRANCH ${AMANZI_GIT_STATUS})
  endif()

  message(STATUS ">>>> JDM: AMANZI_GIT_BRANCH = ${AMANZI_GIT_BRANCH}")

  # Extract the lastest tag of the form amanzi-*

  # Get the hash of the current version
  set(GIT_ARGS rev-parse --short HEAD)
  execute_process(COMMAND  ${GIT_EXECUTABLE} ${GIT_ARGS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  RESULT_VARIABLE err_occurred
                  OUTPUT_VARIABLE AMANZI_GIT_GLOBAL_HASH
                  ERROR_VARIABLE err
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)
  if(err_occurred)
    message(WARNING "Error executing git:\n ${cmd}\n${err}")
    set(cmd_output cmd_output-NOTFOUND)
    exit()
  endif()

  message(STATUS ">>>> JDM: AMANZI_GIT_GLOBAL_HASH: ${AMANZI_GIT_GLOBAL_HASH}")

  # Ensure repository has the latest tags
  set(GIT_ARGS fetch --no-recurse-submodules --tags)
  execute_process(COMMAND  ${GIT_EXECUTABLE} ${GIT_ARGS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  RESULT_VARIABLE err_occurred
                  OUTPUT_VARIABLE cmd_output
                  ERROR_VARIABLE err
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)
  if(err_occurred)
    message(WARNING "Error executing git:\n ${cmd}\n${cmd_output}\n${err}")
    set(cmd_output cmd_output-NOTFOUND)
    exit()
  endif()

  # Get the latest amanzi-* version number tag
  set(GIT_ARGS tag -l amanzi-*)
  execute_process(COMMAND  ${GIT_EXECUTABLE} ${GIT_ARGS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  RESULT_VARIABLE err_occurred
                  OUTPUT_VARIABLE AMANZI_GIT_LATEST_TAG
                  ERROR_VARIABLE err
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)

  # Put the tags in a list
  STRING(REPLACE "\n" ";" AMANZI_GIT_LATEST_TAG_LIST ${AMANZI_GIT_LATEST_TAG})
  # Extract the lastest tag of the form amanzi-*
  IF ( ${AMANZI_GIT_BRANCH} MATCHES "master" )
    FOREACH(atag ${AMANZI_GIT_LATEST_TAG_LIST})
      IF ( ${atag} MATCHES "^amanzi-.*-dev" )
        set ( AMANZI_GIT_LATEST_TAG ${atag} )
      ENDIF()
    ENDFOREACH()
  ELSE()
    FOREACH(atag ${AMANZI_GIT_LATEST_TAG_LIST})
      IF ( ${atag} MATCHES "^amanzi-[0-9]\\.[0-9][0-9]?\\.[0-9]" )
        set ( AMANZI_GIT_LATEST_TAG ${atag} )
      ENDIF()
    ENDFOREACH()
  ENDIF()

  message(STATUS ">>>> JDM: GIT_EXEC        = ${GIT_EXECUTABLE}")
  message(STATUS ">>>> JDM: GIT_ARGS        = ${GIT_ARGS}")
  message(STATUS ">>>> JDM: RESULT_VARIABLE = ${err_occurred}")
  message(STATUS ">>>> JDM: AMANZI_GIT_LATEST_TAG = ${AMANZI_GIT_LATEST_TAG}")

  #STRING(REGEX REPLACE "amanzi-" "" AMANZI_GIT_LATEST_TAG_VER ${AMANZI_GIT_LATEST_TAG})
  #STRING(REGEX REPLACE "\\..*" "" AMANZI_GIT_LATEST_TAG_MAJOR ${AMANZI_GIT_LATEST_TAG_VER})
  #STRING(REGEX MATCH "\\.[0-9][0-9]?[\\.,-]" AMANZI_GIT_LATEST_TAG_MINOR ${AMANZI_GIT_LATEST_TAG_VER})
  #STRING(REGEX REPLACE "[\\.,-]" "" AMANZI_GIT_LATEST_TAG_MINOR ${AMANZI_GIT_LATEST_TAG_MINOR} )

  #set(AMANZI_VERSION_MAJOR ${AMANZI_GIT_LATEST_TAG_MAJOR})
  #set(AMANZI_VERSION_MINOR ${AMANZI_GIT_LATEST_TAG_MINOR})

  #message(STATUS ">>>> JDM: AMANZI_GIT_LATEST_TAG_MAJOR = ${AMANZI_GIT_LATEST_TAG_MAJOR}")
  #message(STATUS ">>>> JDM: AMANZI_GIT_LATEST_TAG_MINOR = ${AMANZI_GIT_LATEST_TAG_MINOR}")

  #
  # Amanzi version
  #
  set(AMANZI_VERSION ${AMANZI_GIT_LATEST_TAG_VER}_${AMANZI_GIT_GLOBAL_HASH})

  #STRING(REGEX REPLACE ".*\\.[0-9][0-9]?[\\.,-]" "" AMANZI_VERSION_PATCH ${AMANZI_VERSION})
  #STRING(REGEX REPLACE ".*_" "" AMANZI_VERSION_HASH ${AMANZI_VERSION_PATCH})
  #STRING(REGEX REPLACE "_.*" "" AMANZI_VERSION_PATCH ${AMANZI_VERSION_PATCH})

else()

  message(STATUS "  >>>>>>>> Using static version information to create amanzi_version.hh")
  if ( NOT GIT_FOUND )
    message(STATUS "    >>>>>> Could not locate Git executable.")
  endif()
  if ( NOT EXISTS ${CMAKE_SOURCE_DIR}/.git/ )
    message(STATUS "    >>>>>> Release or snapshot, no .git directory found.")
  endif()

  #
  # For releases and snapshots, set static information before creating the tarball.
  #
  set(AMANZI_GIT_BRANCH master )
  set(AMANZI_GIT_GLOBAL_HASH )

  set(AMANZI_VERSION_MAJOR 1)
  set(AMANZI_VERSION_MINOR 0)
  set(AMANZI_VERSION_PATCH 0)
  set(AMANZI_VERSION_HASH ${AMANZI_GIT_GLOBAL_HASH})

  #
  # Amanzi version
  #
  set(AMANZI_VERSION ${AMANZI_VERSION_MAJOR}.${AMANZI_VERSION_MINOR}-${AMANZI_VERSION_PATCH}_${AMANZI_VERSION_HASH})

endif()

# Write the version header file
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
