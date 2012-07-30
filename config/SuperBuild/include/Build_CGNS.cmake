#  -*- mode: cmake -*-

#
# Build TPL:  CGNS 
#  

# --- Define all the directories and common external project flags
define_external_project_args(CGNS
                             TARGET cgns
                             BUILD_IN_SOURCE)

# --- Define the configure command
set(CGNS_cmake_configure ${CGNS_prefix_dir}/cgns-configure-step.cmake)
set(CGNS_sh_configure    ${CGNS_prefix_dir}/cgns-configure-step.sh)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/cgns-configure-step.cmake.in
               ${CGNS_cmake_configure}
               @ONLY)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/cgns-configure-step.sh.in
               ${CGNS_sh_configure}
               @ONLY)
set(CGNS_CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${CGNS_cmake_configure})



# --- Define the install command
set(CGNS_cmake_install ${CGNS_prefix_dir}/cgns-install-step.cmake)
set(CGNS_sh_install    ${CGNS_prefix_dir}/cgns-install-step.sh)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/cgns-install-step.cmake.in
               ${CGNS_cmake_install}
               @ONLY)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/cgns-install-step.sh.in
               ${CGNS_sh_install}
               @ONLY)
set(CGNS_INSTALL_COMMAND ${CMAKE_COMMAND} -P ${CGNS_cmake_install})

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
                    INSTALL_DIR      ${CGNS_prefix_dir}        # Install directory
                    INSTALL_COMMAND  ${CGNS_INSTALL_COMMAND}   # Install command
                    # -- Output control
                    ${CGNS_logging_args})
