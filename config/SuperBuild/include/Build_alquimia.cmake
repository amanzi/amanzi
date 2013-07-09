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
                    INSTALL_COMMAND   ""
#                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}        # Install directory
                    # -- Output control
                    ${ALQUIMIA_logging_args})

include(BuildLibraryName)
build_library_name(alquimia ALQUIMIA_LIBRARIES APPEND_PATH ${TPL_INSTALL_PREFIX}/lib)
set(ALQUIMIA_INCLUDE_DIRS ${TPL_INSTALL_PREFIX}/include)

