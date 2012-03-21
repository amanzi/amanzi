#  -*- mode: cmake -*-

#
# Build TPL: CURL 
#

# --- Define all the directories and common external project flags
define_external_project_args(CURL
                             TARGET curl
                             DEPENDS ZLIB)

# --- Define the configuration parameters

# Search for OpenSSL and update environment variables 
# so that the CURL build system detects OpenSSL
set(ssl_pkg_config)
find_package(OpenSSL)
if ( OPENSSL_FOUND )
  list(APPEND curl_cflags_list -I${OPENSSL_INCLUDE_DIR})
  foreach(lib ${OPENSSL_LIBRARIES} )
    get_filename_component(_lib_name ${lib} NAME_WE)
    get_filename_component(_lib_ext ${lib}  EXT)
    get_filename_component(_lib_path ${lib} PATH)

    set(_lib_prefix)
    if ( ${_lib_ext} STREQUAL ${CMAKE_SHARED_LIBRARY_SUFFIX} )
      set(_lib_prefix ${CMAKE_SHARED_LIBRARY_PREFIX})
    elseif( ${_lib_ext} STREQUAL ${CMAKE_STATIC_LIBRARY_SUFFIX} )
      set(_lib_prefix ${CMAKE_STATIC_LIBRARY_PREFIX})
    endif()

    if ( ${_lib_name} STREQUAL ${_lib_prefix}ssl ) 
      set(ssl_pkg_config PKG_CONFIG_PATH=${_lib_path})
    endif()

  endforeach()

else()

  message(WARNING "Failed to locate OpenSSL. Curl build may fail."
                  "If it does fail, rerun cmake configuration with"
                  "\n-D OPENSSL_ROOT_DIR:FILEPATH=/OpenSSL/install/prefix\n"
                  "to define the OpenSSL installation prefix.")
endif()

# Now build the CPPFLAGS string WHITESPACE is needed here!
include(BuildWhitespaceString)
build_whitespace_string(curl_cflags -I${TPL_INSTALL_PREFIX}/include ${Amanzi_COMMON_CFLAGS})

# --- Add external project build and tie to the CURL build target
ExternalProject_Add(${CURL_BUILD_TARGET}
                    DEPENDS   ${CURL_PACKAGE_DEPENDS}             # Package dependency target
                    TMP_DIR   ${CURL_tmp_dir}                     # Temporary files directory
                    STAMP_DIR ${CURL_stamp_dir}                   # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}              # Download directory
                    URL          ${CURL_URL}                      # URL may be a web site OR a local file
                    URL_MD5      ${CURL_MD5_SUM}                  # md5sum of the archive file
                    # -- Configure
                    SOURCE_DIR       ${CURL_source_dir}           # Source directory
                    CONFIGURE_COMMAND
                                     <SOURCE_DIR>/configure
                                         --prefix=<INSTALL_DIR>
                                         --enable-static
                                         --disable-shared
                                         --with-zlib=${TPL_INSTALL_PREFIX}/lib
                                         CC=${CMAKE_C_COMPILER}
                                         CFLAGS=${curl_cflags}
                                         LDFLAGS=-L${TPL_INSTALL_PREFIX}/lib
                                         ${ssl_pkg_config}
                    # -- Build
                    BINARY_DIR        ${CURL_build_dir}           # Build directory 
                    BUILD_COMMAND     $(MAKE)                     # $(MAKE) enables parallel builds through make
                    BUILD_IN_SOURCE   ${CURL_BUILD_IN_SOURCE}     # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}        # Install directory
                    # -- Output control
                    ${CURL_logging_args})

# --- Define the CURL executable  
global_set(CURL_EXECUTABLE ${TPL_INSTALL_PREFIX}/bin/curl)
