#  -*- mode: cmake -*-

#
# Build TPL:  ZLIB 
#   

# --- Define all the directories and common external project flags
define_external_project_args(ZLIB
                             TARGET zlib
                             BUILD_IN_SOURCE)


# --- Define the CMake configure parameters
# Note:
#      CMAKE_CACHE_ARGS requires -DVAR:<TYPE>=VALUE syntax
#      CMAKE_ARGS -DVAR=VALUE OK
# NO WHITESPACE between -D and VAR. Parser blows up otherwise.
set(ZLIB_CMAKE_CACHE_ARGS
                  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
                  -DCMAKE_INSTALL_PREFIX:STRING=<INSTALL_DIR>
                  -DBUILD_SHARED_LIBS:BOOL=FALSE)

# --- Add external project build and tie to the ZLIB build target
ExternalProject_Add(${ZLIB_BUILD_TARGET}
                    DEPENDS   ${ZLIB_PACKAGE_DEPENDS}             # Package dependency target
                    TMP_DIR   ${ZLIB_tmp_dir}                     # Temporary files directory
                    STAMP_DIR ${ZLIB_stamp_dir}                   # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}              # Download directory
                    URL          ${ZLIB_URL}                      # URL may be a web site OR a local file
                    URL_MD5      ${ZLIB_MD5_SUM}                  # md5sum of the archive file
                    # -- Configure
                    SOURCE_DIR       ${ZLIB_source_dir}               # Source directory
                    CMAKE_CACHE_ARGS ${ZLIB_CMAKE_CACHE_ARGS}         # CMAKE_CACHE_ARGS or CMAKE_ARGS => CMake configure
                                     ${Amanzi_CMAKE_C_COMPILER_ARGS}  # Ensure uniform build
                    # -- Build
                    BINARY_DIR        ${ZLIB_build_dir}           # Build directory 
                    BUILD_COMMAND     $(MAKE)                     # $(MAKE) enables parallel builds through make
                    BUILD_IN_SOURCE   ${ZLIB_BUILD_IN_SOURCE}     # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}        # Install directory
                    # -- Output control
                    ${ZLIB_logging_args})
