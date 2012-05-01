#  -*- mode: cmake -*-

#
# Build TPL:  CGNS 
#  

# --- Define all the directories and common external project flags
define_external_project_args(CGNS
                             TARGET cgns
			     DEPENDS ZLIB
                             BUILD_IN_SOURCE)

# --- Define the configure command

# Build the configure script
set(CGNS_sh_configure ${CGNS_prefix_dir}/cgns-configure-step.sh)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/cgns-configure-step.sh.in
               ${CGNS_sh_configure}
	       @ONLY)

# Configure the CMake command file
set(CGNS_cmake_configure ${CGNS_prefix_dir}/cgns-configure-step.cmake)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/cgns-configure-step.cmake.in
               ${CGNS_cmake_configure}
	       @ONLY)
set(CGNS_CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${CGNS_cmake_configure})	

# --- Add external project build and tie to the CGNS build target
ExternalProject_Add(${CGNS_BUILD_TARGET}
                    DEPENDS   ${CGNS_PACKAGE_DEPENDS}             # Package dependency target
                    TMP_DIR   ${CGNS_tmp_dir}                     # Temporary files directory
                    STAMP_DIR ${CGNS_stamp_dir}                   # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}              # Download directory
                    URL          ${CGNS_URL}                      # URL may be a web site OR a local file
                    URL_MD5      ${CGNS_MD5_SUM}                  # md5sum of the archive file
                    # -- Configure
                    SOURCE_DIR        ${CGNS_source_dir}               # Source directory
                    CONFIGURE_COMMAND ${CGNS_CONFIGURE_COMMAND}
                    # -- Build
                    BINARY_DIR        ${CGNS_build_dir}           # Build directory 
                    BUILD_COMMAND     $(MAKE)                     # $(MAKE) enables parallel builds through make
                    BUILD_IN_SOURCE   ${CGNS_BUILD_IN_SOURCE}     # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}        # Install directory
                    # -- Output control
                    ${CGNS_logging_args})

# Add the external project and tie to zlib target
#ExternalProject_Add(${CGNS_target}
#    DEPENDS curl
#    ${CGNS_ep_directory_args}
#    ${CGNS_url_args}
#    CONFIGURE_COMMAND
#                    <SOURCE_DIR>/configure 
#                              --prefix=<INSTALL_DIR>
#                              --enable-lfs
#                              --enable-64bit
#                              --with-fortran=no
#                              FC=false
#    BUILD_IN_SOURCE TRUE      
#    ${CGNS_logging_opts}               
#)
