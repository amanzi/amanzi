#  -*- mode: cmake -*-

#
# Build TPL: ASCEMIO 
#  

# --- Define all the directories and common external project flags
define_external_project_args(ASCEMIO 
                             TARGET ascemio
                             BUILD_IN_SOURCE
                             DEPENDS HDF5)


# -- Define the build command

# Build the build script
set(ASCEMIO_sh_build ${ASCEMIO_prefix_dir}/ascemio-build-step.sh)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/ascemio-build-step.sh.in
               ${ASCEMIO_sh_build}
               @ONLY)

# Configure the CMake command file
set(ASCEMIO_cmake_build ${ASCEMIO_prefix_dir}/ascemio-build-step.cmake)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/ascemio-build-step.cmake.in
               ${ASCEMIO_cmake_build}
               @ONLY)
set(ASCEMIO_BUILD_COMMAND ${CMAKE_COMMAND} -P ${ASCEMIO_cmake_build})

# --- Define the install command

# Build the install script
set(ASCEMIO_sh_install ${ASCEMIO_prefix_dir}/ascemio-install-step.sh)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/ascemio-install-step.sh.in
               ${ASCEMIO_sh_install}
               @ONLY)

# Configure the CMake command file
set(ASCEMIO_cmake_install ${ASCEMIO_prefix_dir}/ascemio-install-step.cmake)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/ascemio-install-step.cmake.in
               ${ASCEMIO_cmake_install}
               @ONLY)
set(ASCEMIO_INSTALL_COMMAND ${CMAKE_COMMAND} -P ${ASCEMIO_cmake_install})

# --- Add external project build and tie to the ZLIB build target
ExternalProject_Add(${ASCEMIO_BUILD_TARGET}
                    DEPENDS   ${ASCEMIO_PACKAGE_DEPENDS}             # Package dependency target
                    TMP_DIR   ${ASCEMIO_tmp_dir}                     # Temporary files directory
                    STAMP_DIR ${ASCEMIO_stamp_dir}                   # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}                 # Download directory
                    URL          ${ASCEMIO_URL}                      # URL may be a web site OR a local file
                    URL_MD5      ${ASCEMIO_MD5_SUM}                  # md5sum of the archive file
                    # -- Configure
                    SOURCE_DIR       ${ASCEMIO_source_dir}           # Source directory
                    CONFIGURE_COMMAND ""
                    # -- Build
                    BINARY_DIR        ${ASCEMIO_build_dir}           # Build directory 
                    BUILD_COMMAND     ${ASCEMIO_BUILD_COMMAND}       # Build command
                    BUILD_IN_SOURCE   ${ASCEMIO_BUILD_IN_SOURCE}     # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}           # Install directory
                    INSTALL_COMMAND  ${ASCEMIO_INSTALL_COMMAND}      # Install command
                    # -- Output control
                    ${ASCEMIO_logging_args})
