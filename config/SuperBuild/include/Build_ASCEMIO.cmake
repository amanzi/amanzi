#  -*- mode: cmake -*-

#
# Build TPL: ASCEMIO 
#  

# --- Define all the directories and common external project flags
define_external_project_args(ASCEMIO 
                             TARGET ascemio
                             BUILD_IN_SOURCE
			     DEPENDS HDF5)


message("IN HERE")			   
#DEBUG
#DEBUG# ########################################################################### #
#DEBUG# Build the build command
#DEBUG# ########################################################################### #
#DEBUG
#DEBUG# Build the build script
#DEBUGset(ASCEMIO_sh_build ${ASCEMIO_prefix_dir}/ascemio-build-step.sh)
#DEBUGconfigure_file(${SuperBuild_BUILD_FILES_DIR}/ascemio-build-step.sh.in
#DEBUG               ${ASCEMIO_sh_build}
#DEBUG	       @ONLY)
#DEBUG
#DEBUG# Configure the CMake command file
#DEBUGset(ASCEMIO_cmake_build ${ASCEMIO_prefix_dir}/ascemio-build-step.cmake)
#DEBUGconfigure_file(${SuperBuild_BUILD_FILES_DIR}/ascemio-build-step.cmake.in
#DEBUG               ${ASCEMIO_cmake_build}
#DEBUG	       @ONLY)
#DEBUGset(ASCEMIO_BUILD_COMMAND ${CMAKE_COMMAND} -P ${ASCEMIO_cmake_build})	
#DEBUG
#DEBUG# ########################################################################### #
#DEBUG# Build the install command
#DEBUG# ########################################################################### #
#DEBUG
#DEBUG# Build the install script
#DEBUGset(ASCEMIO_sh_install ${ASCEMIO_prefix_dir}/ascemio-install-step.sh)
#DEBUGconfigure_file(${SuperBuild_BUILD_FILES_DIR}/ascemio-install-step.sh.in
#DEBUG               ${ASCEMIO_sh_install}
#DEBUG	       @ONLY)
#DEBUG
#DEBUG# Configure the CMake command file
#DEBUGset(ASCEMIO_cmake_install ${ASCEMIO_prefix_dir}/ascemio-install-step.cmake)
#DEBUGconfigure_file(${SuperBuild_BUILD_FILES_DIR}/ascemio-install-step.cmake.in
#DEBUG               ${ASCEMIO_cmake_install}
#DEBUG	       @ONLY)
#DEBUGset(ASCEMIO_INSTALL_COMMAND ${CMAKE_COMMAND} -P ${ASCEMIO_cmake_install})	
#DEBUG
#DEBUG# ########################################################################### #
#DEBUG# Add External Project
#DEBUG# ########################################################################### #
#DEBUGExternalProject_Add(${ASCEMIO_target}
#DEBUG    DEPENDS hdf5
#DEBUG    ${ASCEMIO_ep_directory_args}
#DEBUG    ${ASCEMIO_url_args}
#DEBUG    CONFIGURE_COMMAND ""
#DEBUG    BUILD_COMMAND ${ASCEMIO_BUILD_COMMAND}
#DEBUG    INSTALL_COMMAND ${ASCEMIO_INSTALL_COMMAND}
#DEBUG    ${ASCEMIO_logging_opts}
#DEBUG)
