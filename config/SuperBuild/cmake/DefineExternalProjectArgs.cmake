# ############################################################################ #
#                                                                              #
#  DEFINE_EXTERNAL_PROJECT_ARGS(<PACK_NAME>                                    # 
#                                [ TARGET target-name ]                        #
#                                [ BUILD_IN_SOURCE ]                           #
#                                [ DEPENDS pack1 pack2 pack3 ]                 #
#                              )                                               #
#                                                                              #
#  A macro that defines common arguments for the AddExternalProject function   #
#  This macro provides an organized build structure.                           #
#                                                                              #
# ############################################################################ #
include(CMakeParseArguments)
include(SetMacros)
macro(DEFINE_EXTERNAL_PROJECT_ARGS prefix)
 
  # --- Parse the arguments
  set(_flags "BUILD_IN_SOURCE")
  set(_oneValue "TARGET")
  set(_multiValue "DEPENDS")
  cmake_parse_arguments(PARSE "${_flags}" "${_oneValue}" "${_multiValue}" ${ARGN})


  # --- Define the build target name
  if ( NOT PARSE_TARGET )
    string(TOLOWER "${prefix}" _target_name)
  else()
    set(_target_name ${PARSE_TARGET})
  endif()
  global_set(${prefix}_BUILD_TARGET ${_target_name})


  # --  Define the directories for download, build, install and timestamps

  # Will use lower case for directory names
  string(TOLOWER "${${prefix}_BUILD_TARGET}" target_lc) 


  set(${prefix}_prefix_dir ${SuperBuild_BINARY_DIR}/${target_lc})
  set(${prefix}_source_dir ${SuperBuild_BINARY_DIR}/${target_lc}/${target_lc}-${${prefix}_VERSION}-source)
  set(${prefix}_stamp_dir  ${SuperBuild_BINARY_DIR}/${target_lc}/${target_lc}-timestamps)
  set(${prefix}_tmp_dir    ${SuperBuild_BINARY_DIR}/${target_lc}/tmp)

  # Default is to build out of source, but some packages can not do that
  if ( NOT PARSE_BUILD_IN_SOURCE ) 
    set(${prefix}_build_dir  ${SuperBuild_BINARY_DIR}/${target_lc}/${target_lc}-${${prefix}_VERSION}-build)
  else()  
    set(${prefix}_build_dir "")
  endif()

  # Download from the web unless DISABLE_EXTERNAL_DOWNLOADS is TRUE
  if ( DISABLE_EXTERNAL_DOWNLOAD )
    message(STATUS ">>>> JDM: Project = ${prefix}")
    message(STATUS ">>>> JDM: GIT REPOSITORY = ${${prefix}_GIT_REPOSITORY}")
    if ( DEFINED ${prefix}_GIT_REPOSITORY )
      if ( ${${prefix}_GIT_REPOSITORY} MATCHES "https://bitbucket.*" )
        message(STATUS ">>>> JDM: This is a bitbucket repository.")
        STRING(REGEX REPLACE "/src/master" "" GIT_REPOSITORY_TEMP ${${prefix}_GIT_REPOSITORY})
        message(STATUS ">>>> JDM: GIT_REPOSITORY_TEMP = ${GIT_REPOSITORY_TEMP}")
        STRING(REGEX REPLACE ".*\/" "" GIT_REPOSITORY_LOCAL_DIR ${GIT_REPOSITORY_TEMP})
      else()
        STRING(REGEX REPLACE ".*\/" "" GIT_REPOSITORY_LOCAL_DIR ${${prefix}_GIT_REPOSITORY})
      endif()
      message(STATUS ">>>> JDM: REGEX = ${GIT_REPOSITORY_LOCAL_DIR}")
      set(${prefix}_URL ${TPL_DOWNLOAD_DIR}/${GIT_REPOSITORY_LOCAL_DIR})
    else()
	  set(${prefix}_URL ${TPL_DOWNLOAD_DIR}/${${prefix}_SAVEAS_FILE})
    endif()
    
    if ( NOT EXISTS "${${prefix}_URL}" )
      message(FATAL_ERROR "You have disabled external downloads (-DDISABLE_EXTERNAL_DOWNLOAD:BOOL=TRUE),"
	                  "however ${${prefix}_URL} does not exist locally on this system.")
    endif()

  else()
    set(${prefix}_URL ${${prefix}_URL_STRING}/${${prefix}_ARCHIVE_FILE})
  endif()

  # --- Set additional arguments

  # Log all steps, this keeps the STDOUT/STDERR tidy.
  set(${prefix}_logging_args
                          LOG_DOWNLOAD  1 
                          LOG_UPDATE    1 
                          LOG_CONFIGURE 1
                          LOG_BUILD     1
                          LOG_TEST      1
                          LOG_INSTALL   1)

  # Define the package dependencies 			
  set(${prefix}_PACKAGE_DEPENDS)
  foreach( _pack ${PARSE_DEPENDS})
    set(_pack_target "${${_pack}_BUILD_TARGET}")
    if ( NOT TARGET ${_pack_target} )
      message(FATAL_ERROR "Package ${prefix} requires ${_pack}, "
	                  "however the build configuration for ${_pack} has not been defined.")
    endif()
    list(APPEND ${prefix}_PACKAGE_DEPENDS "${_pack_target}")
  endforeach()

  # Set the build in source flag
  set(${prefix}_BUILD_IN_SOURCE ${PARSE_BUILD_IN_SOURCE})

endmacro(DEFINE_EXTERNAL_PROJECT_ARGS)
