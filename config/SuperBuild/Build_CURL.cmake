#  -*- mode: cmake -*-

#
# Build CURL 
#    

include(ExternalProject)
include(BuildLibraryName)

include(TPLVersions)


# Define source, build and install directories
set(CURL_source_dir "${CMAKE_BINARY_DIR}/external-projects/curl/curl-${CURL_VERSION}-source")
set(CURL_binary_dir "${CMAKE_BINARY_DIR}/external-projects/curl/curl-${CURL_VERSION}-build")
set(CURL_install_dir "${CMAKE_BINARY_DIR}/external-projects/curl")

# Make target: build curl with 'make curl' command
set(CURL_target curl)

# Add the external project and tie to curl target
ExternalProject_Add(${CURL_target}
    SOURCE_DIR ${CURL_source_dir}
    BINARY_DIR ${CURL_binary_dir}
    INSTALL_DIR ${CURL_install_dir}
    URL ${CURL_URL_STRING}/${CURL_ARCHIVE_FILE}
    URL_MD5 ${CURL_MD5_SUM}
    CONFIGURE_COMMAND
                   <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --enable-static
)
                                      
# Define variables needed by other packages
# These should match the final output from FindCURL.cmake

# Include directories
set(CURL_INCLUDE_DIR "${CURL_install_dir}/include")
set(CURL_INCLUDE_DIRS "${CURL_install_dir}/include")

# Libraries
build_library_name(curl CURL_LIBRARY_FILENAME STATIC)
set(CURL_LIBRARY ${CURL_install_dir}/lib/${CURL_LIBRARY_FILENAME})
set(CURL_LIBRARIES "${CURL_LIBRARY}")
