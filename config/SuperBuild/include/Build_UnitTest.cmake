#  -*- mode: cmake -*-

#
# Build TPL: UnitTest
# 

# --- Define all the directories and common external project flags
define_external_project_args(UnitTest
                             TARGET unittest
			     BUILD_IN_SOURCE)

# --- Define patch command

# Need Perl to patch
find_package(Perl)
if ( NOT PERL_FOUND )
  message(FATAL_ERROR "Failed to locate perl. "
                      "Can not patch UnitTest without PERL")
endif()

# Build the patch script
set(UnitTest_sh_patch ${UnitTest_prefix_dir}/unittest-patch-step.sh)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/unittest-patch-step.sh.in
               ${UnitTest_sh_patch}
               @ONLY)

# --- Define the install command

# Build the install script
set(UnitTest_sh_install ${UnitTest_prefix_dir}/unittest-install-step.sh)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/unittest-install-step.sh.in
               ${UnitTest_sh_install}
               @ONLY)

	     
# --- Add external project build and tie to the ZLIB build target
ExternalProject_add(${UnitTest_BUILD_TARGET}
                    DEPENDS   ${UnitTest_PACKAGE_DEPENDS}             # Package dependency target
                    TMP_DIR   ${UnitTest_tmp_dir}                     # Temporary files directory
                    STAMP_DIR ${UnitTest_stamp_dir}                   # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}                  # Download directory
                    URL          ${UnitTest_URL}                      # URL may be a web site OR a local file
                    URL_MD5      ${UnitTest_MD5_SUM}                  # md5sum of the archive file
                    # -- Patch
		    PATCH_COMMAND sh ${UnitTest_sh_patch}             # Run the patch script
		    # -- Configure
		    CONFIGURE_COMMAND   ""                            # No configure step
		    SOURCE_DIR          ${UnitTest_source_dir}        # Defining forces CMake to mkdir SOURCE_DIR
		    # -- Build
		    BUILD_COMMAND       $(MAKE)                       # Run make in build directory $(MAKE) enables parallel build
		    BINARY_DIR          ${UnitTest_build_dir}         # Define the build directory
		    BUILD_IN_SOURCE     ${UnitTest_BUILD_IN_SOURCE}   # Flag in/out source build
                    # -- Install
                    INSTALL_DIR         ${TPL_INSTALL_PREFIX}        # Install directory
		    INSTALL_COMMAND     sh ${UnitTest_sh_install}    # Run the install script
                    # -- Output control
                    ${UnitTest_logging_args})
