#  -*- mode: cmake -*-

#
# Build TPL: Boost 
#    
define_external_project(Boost TARGET boost)

# ########################################################################### #
# General build definitions
# ########################################################################### #

# We only build what we need, this is NOT a full Boost install
set(Boost_projects "system,filesystem,program_options,regex")

# Boost download URL is not a standard <url>/<filename>.tar.gz pattern
# User must download this file in a local directory, until I find a fix
#if(NOT TPL_Boost_DOWNLOAD_DIR)
#  set(TPL_Boost_DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR})
#endif()
#set(Boost_url ${TPL_Boost_DOWNLOAD_DIR}/${Boost_ARCHIVE_FILE})

#if ( NOT EXISTS ${Boost_url} )
#  message(FATAL_ERROR "Download of Boost not supported at this time. "
#                      "Please define the location of ${Boost_ARCHIVE_FILE} "
#		      "with -D TPL_Boost_DOWNLOAD_DIR:FILEPATH=<some path>.")
#endif()		    

# ########################################################################### #
# Build the configure command
# ########################################################################### #

# Determine toolset type
set(Boost_toolset)
string(TOLOWER ${CMAKE_C_COMPILER_ID} compiler_id_lc)
if (compiler_id_lc)
  if (APPLE)
    if ( ${compiler_id_lc} STREQUAL "intel" )
      set(Boost_toolset --with-toolset=intel-darwin)
    else()  
      set(Boost_toolset --with-toolset=darwin)
    endif()  
  elseif(UNIX)
    if ( ${compiler_id_lc} STREQUAL "gnu" )
        set(Boost_toolset --with-toolset=gcc)
    elseif(${compiler_id_lc} STREQUAL "intel")
        set(Boost_toolset --with-toolset=intel-linux)
    elseif(${compiler_id_lc} STREQUAL "pgi")
        set(Boost_toolset --with-toolset=pgi)
    elseif(${compiler_id_lc} STREQUAL "pathscale")
        set(Boost_toolset --with-toolset=pathscale)
    endif()
  endif()
endif()

configure_file(${SuperBuild_BUILD_FILES_DIR}/boost-configure-step.cmake.in
               ${Boost_prefix_dir}/boost-configure-step.cmake
	       @ONLY)
set(Boost_CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${Boost_prefix_dir}/boost-configure-step.cmake)

# ########################################################################### #
# Build the build command
# ########################################################################### #

configure_file(${SuperBuild_BUILD_FILES_DIR}/boost-build-step.cmake.in
               ${Boost_prefix_dir}/boost-build-step.cmake
	       @ONLY)

set(Boost_BUILD_COMMAND ${CMAKE_COMMAND} -P ${Boost_prefix_dir}/boost-build-step.cmake)	     

# ########################################################################### #
# Add the Boost build target
# ########################################################################### #

# Make target: build boost with 'make boost' command
set(Boost_target boost)

# Add the external project and tie to boost target
ExternalProject_Add(${Boost_target}
    ${Boost_ep_directory_args}
    ${Boost_url_args}
    CONFIGURE_COMMAND ${Boost_CONFIGURE_COMMAND}
    BUILD_COMMAND ${Boost_BUILD_COMMAND}
    INSTALL_COMMAND ""
    ${Boost_logging_opts}
)

# ########################################################################### #
# Define the variables needed to build against this installation
# ########################################################################### #

# Flag find modules
#global_set(Boost_FOUND TRUE)

# Include directories
#global_set(Boost_INCLUDE_DIR "${Boost_install_dir}/include")
#global_set(Boost_INCLUDE_DIRS "${Boost_INCLUDE_DIR}")

# Project Libraries
#string(REGEX REPLACE "," ";" Boost_projects_list ${Boost_projects})
#set(libraries)
#foreach ( proj ${Boost_projects_list} )
#  set(lib boost_${proj})
#  build_library_name(${lib} library)
#  string(TOUPPER ${proj} proj_uc)
#  global_set(Boost_${proj_uc}_LIBRARY ${Boost_install_dir}/lib/${library})
#  global_set(Boost_${proj_uc}_FOUND TRUE)
#  add_library(${lib} UNKNOWN IMPORTED)
#  set_property(TARGET ${lib} 
#               PROPERTY IMPORTED_LOCATION ${Boost_${proj_uc}_LIBRARY})
#  add_dependencies(${lib} ${Boost_target})	     

#  #list(APPEND libraries ${Boost_${proj_uc}_LIBRARY})
#  list(APPEND libraries ${lib})
#endforeach() 
#global_set(Boost_LIBRARIES ${libraries})
#global_set(Boost_LIBRARIES boost_system)
#if ( NOT TARGET boost_system )
#  message(FATAL_ERROR "boost_system is not a target")
#else ()
#    message(STATUS "boost_system is a target")
#endif()    




