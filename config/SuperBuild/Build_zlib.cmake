#  -*- mode: cmake -*-

#
# Build ZLIB 
#    
include(ExternalProject)
include(BuildLibraryName)

include(TPLVersions)


# Define source, build and install directories
# Out of source builds broken in zlib. No binary directory defined.
set(ZLIB_source_dir "${CMAKE_BINARY_DIR}/external-projects/zlib/zlib-${ZLIB_VERSION}-source")
set(ZLIB_install_dir "${CMAKE_BINARY_DIR}/external-projects/zlib")

# Make target: build zlib with 'make zlib' command
set(ZLIB_target zlib)

# Add the external project and tie to zlib target
ExternalProject_Add(${ZLIB_target}
    SOURCE_DIR ${ZLIB_source_dir}
    BINARY_DIR ${ZLIB_binary_dir}
    INSTALL_DIR ${ZLIB_install_dir}
    URL ${ZLIB_URL_STRING}/${ZLIB_ARCHIVE_FILE}
    URL_MD5 ${ZLIB_MD5_SUM}
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND
                   <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
)
                                      
# Define variables needed by other packages
# These should match the final output from FindZLIB.cmake

# Include directories
set(ZLIB_INCLUDE_DIR "${ZLIB_install_dir}/include")
set(ZLIB_INCLUDE_DIRS "${ZLIB_install_dir}/include")

# Libraries
# Issue on MacOSX when zlib is built with mpicc
# Builds *.so files instead of *.dylib files
# Will use the static library to build other TPLS
build_library_name(z ZLIB_LIBRARY_STATIC_FILENAME STATIC)
build_library_name(z ZLIB_LIBRARY_SHARED_FILENAME SHARED)
set(ZLIB_LIBRARY ${ZLIB_install_dir}/lib/${ZLIB_LIBRARY_STATIC_FILENAME})
set(ZLIB_LIBRARIES "${ZLIB_LIBRARY}")

# Definitions for consistency with FindZLIB.cmake
set(ZLIB_VERSION_STRING ${ZLIB_VERSION})
set(ZLIB_MAJOR_VERSION ${ZLIB_VERSION_MAJOR})
set(ZLIB_MINOR_VERSION ${ZLIB_VERSION_MINOR})
set(ZLIB_PATCH_VERSION ${ZLIB_VERSION_PATCH})
