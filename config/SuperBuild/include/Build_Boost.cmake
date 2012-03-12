#  -*- mode: cmake -*-

#
# Build TPL: Boost 
#

# --- Define all the directories and common external project flags
define_external_project_args(Boost TARGET boost)

# -- Define build definitions

# We only build what we need, this is NOT a full Boost install
set(Boost_projects "system,filesystem,program_options,regex")

# --- Define the configure command

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

configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/boost-configure-step.cmake.in
               ${Boost_prefix_dir}/boost-configure-step.cmake
        @ONLY)
set(Boost_CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${Boost_prefix_dir}/boost-configure-step.cmake)

# --- Define the build command

configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/boost-build-step.cmake.in
               ${Boost_prefix_dir}/boost-build-step.cmake
       @ONLY)

set(Boost_BUILD_COMMAND ${CMAKE_COMMAND} -P ${Boost_prefix_dir}/boost-build-step.cmake)     

# --- Add external project build and tie to the ZLIB build target
ExternalProject_Add(${Boost_BUILD_TARGET}
                    DEPENDS   ${Boost_PACKAGE_DEPENDS}             # Package dependency target
                    TMP_DIR   ${Boost_tmp_dir}                     # Temporary files directory
                    STAMP_DIR ${Boost_stamp_dir}                   # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}              # Download directory
                    URL          ${Boost_URL}                      # URL may be a web site OR a local file
                    URL_MD5      ${Boost_MD5_SUM}                  # md5sum of the archive file
                    # -- Configure
                    SOURCE_DIR       ${Boost_source_dir}           # Source directory
                    CONFIGURE_COMMAND ${Boost_CONFIGURE_COMMAND}
                    # -- Build
                    BINARY_DIR        ${Boost_build_dir}           # Build directory 
                    BUILD_COMMAND     ${Boost_BUILD_COMMAND}       # $(MAKE) enables parallel builds through make
                    BUILD_IN_SOURCE   ${Boost_BUILD_IN_SOURCE}     # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}        # Install directory
                    INSTALL_COMMAND  ""
                    # -- Output control
                    ${Boost_logging_args})
