#  -*- mode: cmake -*-

#
# Build TPL: METIS 
# 

# --- Define all the directories and common external project flags
define_external_project_args(METIS
                             TARGET metis 
			     BUILD_IN_SOURCE)

# --- Define the configure command
find_package(Perl)
if(NOT PERL_FOUND)
  message(FATAL_ERROR "Could not locate Perl. METIS build requires Perl")
endif()

# Build the configure script
set(METIS_sh_configure ${METIS_prefix_dir}/metis-configure-step.sh)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/metis-configure-step.sh.in
               ${METIS_sh_configure}
	       @ONLY)

# Configure the CMake command file
set(METIS_cmake_configure ${METIS_prefix_dir}/metis-configure-step.cmake)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/metis-configure-step.cmake.in
               ${METIS_cmake_configure}
	       @ONLY)
set(METIS_CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${METIS_cmake_configure})	 

# --- Define the build command

# Build the build script
set(METIS_sh_build ${METIS_prefix_dir}/metis-build-step.sh)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/metis-build-step.sh.in
               ${METIS_sh_build}
	       @ONLY)

# Configure the CMake command file
set(METIS_cmake_build ${METIS_prefix_dir}/metis-build-step.cmake)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/metis-build-step.cmake.in
               ${METIS_cmake_build}
	       @ONLY)
set(METIS_BUILD_COMMAND ${CMAKE_COMMAND} -P ${METIS_cmake_build})	

# --- Define the install command

# Build the install script
set(METIS_sh_install ${METIS_prefix_dir}/metis-install-step.sh)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/metis-install-step.sh.in
               ${METIS_sh_install}
	       @ONLY)

# Configure the CMake command file
set(METIS_cmake_install ${METIS_prefix_dir}/metis-install-step.cmake)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/metis-install-step.cmake.in
               ${METIS_cmake_install}
	       @ONLY)
set(METIS_INSTALL_COMMAND ${CMAKE_COMMAND} -P ${METIS_cmake_install})	

# --- Add external project build and tie to the ZLIB build target
set(METIS_GKlib_path ${METIS_source_dir}/GKlib)
ExternalProject_Add(${METIS_BUILD_TARGET}
                    DEPENDS   ${METIS_PACKAGE_DEPENDS}              # Package dependency target
                    TMP_DIR   ${METIS_tmp_dir}                      # Temporary files directory
                    STAMP_DIR ${METIS_stamp_dir}                    # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}                # Download directory
                    URL          ${METIS_URL}                       # URL may be a web site OR a local file
                    URL_MD5      ${METIS_MD5_SUM}                   # md5sum of the archive file
                    # -- Configure
                    SOURCE_DIR        ${METIS_source_dir}           # Source directory
		    CONFIGURE_COMMAND ${METIS_CONFIGURE_COMMAND}
                    # -- Build
                    BINARY_DIR        ${METIS_build_dir}            # Build directory 
                    BUILD_COMMAND     ${METIS_BUILD_COMMAND}        # Build command
                    BUILD_IN_SOURCE   ${METIS_BUILD_IN_SOURCE}      # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}          # Install directory
		    INSTALL_COMMAND  ${METIS_INSTALL_COMMAND}       # Install command
                    # -- Output control
                    ${METIS_logging_args})
# ########################################################################### #
#ExternalProject_Add(${METIS_target}
#    ${METIS_ep_directory_args}
#    ${METIS_url_args}
#    CONFIGURE_COMMAND ${METIS_CONFIGURE_COMMAND}
#    INSTALL_COMMAND ${METIS_INSTALL_COMMAND}
#    ${METIS_logging_opts}
#)

# MSTK needs the full library name
include(BuildLibraryName)
build_library_name(metis METIS_LIBRARY STATIC APPEND_PATH ${TPL_INSTALL_PREFIX}/lib)
