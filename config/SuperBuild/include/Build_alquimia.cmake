#  -*- mode: cmake -*-

#
# Build TPL:  ALQUIMIA 
#   

# Alquimia needs PFlotran.
list(APPEND ALQUIMIA_PACKAGE_DEPENDS ${PFLOTRAN_BUILD_TARGET})

# --- Define all the directories and common external project flags
define_external_project_args(ALQUIMIA
                             TARGET alquimia
                             BUILD_IN_SOURCE)

# --- Define the build command

# Build the build script
set(ALQUIMIA_sh_build ${ALQUIMIA_prefix_dir}/alquimia-build-step.sh)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/alquimia-build-step.sh.in
               ${ALQUIMIA_sh_build}
	       @ONLY)

# Configure the CMake command file
set(ALQUIMIA_cmake_build ${ALQUIMIA_prefix_dir}/alquimia-build-step.cmake)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/alquimia-build-step.cmake.in
               ${ALQUIMIA_cmake_build}
	       @ONLY)
set(ALQUIMIA_CMAKE_COMMAND ${CMAKE_COMMAND} -P ${ALQUIMIA_cmake_build})	

# --- Define the install command

# Build the install script
set(ALQUIMIA_sh_install ${ALQUIMIA_prefix_dir}/alquimia-install-step.sh)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/alquimia-install-step.sh.in
               ${ALQUIMIA_sh_install}
	       @ONLY)

# Configure the CMake command file
set(ALQUIMIA_cmake_install ${ALQUIMIA_prefix_dir}/alquimia-install-step.cmake)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/alquimia-install-step.cmake.in
               ${ALQUIMIA_cmake_install}
	       @ONLY)
set(ALQUIMIA_INSTALL_COMMAND ${CMAKE_COMMAND} -P ${ALQUIMIA_cmake_install})	

# --- Add external project build and tie to the ALQUIMIA build target
ExternalProject_Add(${ALQUIMIA_BUILD_TARGET}
                    DEPENDS   ${ALQUIMIA_PACKAGE_DEPENDS}             # Package dependency target
                    TMP_DIR   ${ALQUIMIA_tmp_dir}                     # Temporary files directory
                    STAMP_DIR ${ALQUIMIAstamp_dir}                   # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}              # Download directory
                    URL          ${ALQUIMIA_URL}                      # URL may be a web site OR a local file
                    URL_MD5      ${ALQUIMIA_MD5_SUM}                  # md5sum of the archive file
                    # -- Configure
                    SOURCE_DIR       ${ALQUIMIA_source_dir}               # Source directory
                    CONFIGURE_COMMAND ""
#                    CMAKE_CACHE_ARGS ${ALQUIMIA_CMAKE_CACHE_ARGS}         # CMAKE_CACHE_ARGS or CMAKE_ARGS => CMake configure
#                                     ${Amanzi_CMAKE_C_COMPILER_ARGS}  # Ensure uniform build
#                                     -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                    # -- Build
#                    BINARY_DIR        ${ALQUIMIA_source_dir}           # Build directory 
                    BUILD_COMMAND     ${ALQUIMIA_CMAKE_COMMAND}                     # $(MAKE) enables parallel builds through make
                    BUILD_IN_SOURCE   1    # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}        # Install directory
            		    INSTALL_COMMAND  ${ALQUIMIA_INSTALL_COMMAND}
                    # -- Output control
                    # -- Output control
                    ${ALQUIMIA_logging_args})

include(BuildLibraryName)
build_library_name(alquimia ALQUIMIA_LIBRARIES APPEND_PATH ${TPL_INSTALL_PREFIX}/lib)
set(ALQUIMIA_INCLUDE_DIRS ${TPL_INSTALL_PREFIX}/include)

