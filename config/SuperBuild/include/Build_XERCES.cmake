#  -*- mode: cmake -*-

#
# Build TPL: XERCES 
#  

# --- Define all the directories and common external project flags
define_external_project_args(XERCES 
                             TARGET xerces
                             BUILD_IN_SOURCE
                             DEPENDS ${MPI_PROJECT} )


# -- Define configure parameters
set(cflag_args "-arch x86_64")
set(cxxflag_args "-arch x86_64")


# --- Add external project build and tie to the ZLIB build target
ExternalProject_Add(${XERCES_BUILD_TARGET}
                    DEPENDS   ${XERCES_PACKAGE_DEPENDS}             # Package dependency target
                    TMP_DIR   ${XERCES_tmp_dir}                     # Temporary files directory
                    STAMP_DIR ${XERCES_stamp_dir}                   # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}                # Download directory
                    URL          ${XERCES_URL}                      # URL may be a web site OR a local file
                    URL_MD5      ${XERCES_MD5_SUM}                  # md5sum of the archive file
                    # -- Configure
                    SOURCE_DIR       ${XERCES_source_dir}           # Source directory
                    CONFIGURE_COMMAND 
		                      <SOURCE_DIR>/configure
				                  --prefix=<INSTALL_DIR>
						  CFLAGS=${cflag_args}
						  CXXFLAGS=${cxxflag_args}
                    # -- Build
                    BINARY_DIR        ${XERCES_build_dir}           # Build directory 
		    BUILD_COMMAND     ${MAKE}                       # Build command
                    BUILD_IN_SOURCE   ${XERCES_BUILD_IN_SOURCE}     # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}          # Install directory
		    INSTALL_COMMAND  ${MAKE}                        # Install command
                    # -- Output control
                    ${XERCES_logging_args})

include(BuildLibraryName)
build_library_name(xerces-c XERCES_LIBRARY APPEND_PATH ${TPL_INSTALL_PREFIX}/lib)
set(XERCES_LIBRARIES ${XERCES_LIBRARY})
set(XERCES_INCLUDE_DIRS ${TPL_INSTALL_PREFIX}/include)
