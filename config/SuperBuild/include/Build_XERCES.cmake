#  -*- mode: cmake -*-
#
# Build TPL: XERCES 
#  
# --- Define all the directories and common external project flags
define_external_project_args(XERCES 
                             TARGET xerces)

# add version version to the autogenerated tpl_versions.h file
amanzi_tpl_version_write(FILENAME ${TPL_VERSIONS_INCLUDE_FILE}
  PREFIX XERCES
  VERSION ${XERCES_VERSION_MAJOR} ${XERCES_VERSION_MINOR} ${XERCES_VERSION_PATCH})

# --- Patch original code
set(XERCES_patch_file xerces-libicu.patch xerces-static-lib-install.patch xerces-cmake-c++17.patch)
set(XERCES_sh_patch ${XERCES_prefix_dir}/xerces-patch-step.sh)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/xerces-patch-step.sh.in
               ${XERCES_sh_patch}
               @ONLY)

# configure the CMake patch step
set(XERCES_cmake_patch ${XERCES_prefix_dir}/xerces-patch-step.cmake)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/xerces-patch-step.cmake.in
               ${XERCES_cmake_patch}
               @ONLY)

set(XERCES_PATCH_COMMAND ${CMAKE_COMMAND} -P ${XERCES_cmake_patch})

# -- Set Xerces configuration options
set(XERCES_CMAKE_ARGS "-DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}")
list(APPEND XERCES_CMAKE_ARGS "-DCMAKE_INSTALL_BINDIR:PATH=${TPL_INSTALL_PREFIX}/bin")
list(APPEND XERCES_CMAKE_ARGS "-DCMAKE_INSTALL_LIBDIR:PATH=${TPL_INSTALL_PREFIX}/lib")
list(APPEND XERCES_CMAKE_ARGS "-DCMAKE_INSTALL_DOCDIR:PATH=${TPL_INSTALL_PREFIX}/doc/xerces")
list(APPEND XERCES_CMAKE_ARGS "-DCMAKE_INSTALL_INCLUDEDIR:PATH=${TPL_INSTALL_PREFIX}/include")
list(APPEND XERCES_CMAKE_ARGS "-Dnetwork:BOOL=FALSE")

# --- Override minimum version
if(CMAKE_MAJOR_VERSION VERSION_EQUAL "4")
  list(APPEND XERCES_CMAKE_ARGS "-DCMAKE_POLICY_VERSION_MINIMUM=3.5")
endif()

# Force OSX to use its CoreServices Framework
if (APPLE)
  list(APPEND XERCES_CMAKE_ARGS "-DCMAKE_HOST_APPLE:BOOL=TRUE")
  list(APPEND XERCES_CMAKE_ARGS "-DXERCES_USE_TRANSCODER_MACOSUNICODECONVERTER:BOOL=TRUE")
endif()

# --- Add external project build and tie to the ZLIB build target
ExternalProject_Add(${XERCES_BUILD_TARGET}
                    DEPENDS   ${XERCES_PACKAGE_DEPENDS}        # Package dependency target
                    TMP_DIR   ${XERCES_tmp_dir}                # Temporary files directory
                    STAMP_DIR ${XERCES_stamp_dir}              # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR  ${TPL_DOWNLOAD_DIR}       
                    URL           ${XERCES_URL}                # URL may be a web site OR a local file
                    URL_MD5       ${XERCES_MD5_SUM}            # md5sum of the archive file
                    PATCH_COMMAND ${XERCES_PATCH_COMMAND} 
                    # -- Configure
                    SOURCE_DIR   ${XERCES_source_dir}          # Source directory
                    CMAKE_ARGS   ${AMANZI_CMAKE_CACHE_ARGS}    # Ensure uniform build
                                 ${XERCES_CMAKE_ARGS}
                                 -DCMAKE_C_FLAGS:STRING=${Amanzi_COMMON_CFLAGS}  # Ensure uniform build
                                 -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                                 -DCMAKE_CXX_FLAGS:STRING=${Amanzi_COMMON_CXXFLAGS}
                                 -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                    # -- Build
                    BINARY_DIR       ${XERCES_build_dir}       # Build directory 
                    BUILD_COMMAND    $(MAKE)
                    BUILD_IN_SOURCE  FALSE                     # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}     # Install directory
                    # -- Output control
                    ${XERCES_logging_args})

include(BuildLibraryName)

build_library_name(xerces-c XERCES_LIBRARY APPEND_PATH ${TPL_INSTALL_PREFIX}/lib)
set(XERCES_LIBRARIES ${XERCES_LIBRARY})
set(XERCES_DIR ${TPL_INSTALL_PREFIX})
set(XERCES_INCLUDE_DIRS ${TPL_INSTALL_PREFIX}/include)
