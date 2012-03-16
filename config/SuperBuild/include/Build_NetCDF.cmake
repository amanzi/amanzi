#  -*- mode: cmake -*-

#
# Build TPL: NetCDF 
# 

# --- Define all the directories and common external project flags
define_external_project_args(NetCDF 
                             TARGET netcdf
                             DEPENDS CURL
                            )

# --- Define the patch command

# Need Perl to patch the files
find_package(Perl)
if (NOT PERL_FOUND)
  message(FATAL_ERROR "Can not locate Perl. Unable to patch and build netCDF")
endif()


# Configure the bash patch script
set(NetCDF_sh_patch ${NetCDF_prefix_dir}/netcdf-patch-step.sh)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/netcdf-patch-step.sh.in
               ${NetCDF_sh_patch}
               @ONLY)

# Configure the CMake command file
set(NetCDF_cmake_patch ${NetCDF_prefix_dir}/netcdf-patch-step.cmake)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/netcdf-patch-step.cmake.in
               ${NetCDF_cmake_patch}
               @ONLY)
set(NetCDF_PATCH_COMMAND ${CMAKE_COMMAND} -P ${NetCDF_cmake_patch})     

# --- Define the configure command

# Default is to build with NetCDF4 which depends on HDF5
option(ENABLE_NetCDF4 "Enable netCDF4 build" TRUE)
set(NetCDF_netcdf4_opts)
if (ENABLE_NetCDF4)

  list(APPEND NetCDF_PACKAGE_DEPENDS ${HDF5_BUILD_TARGET})

  append_set(NetCDF_netcdf4_opts --enable-netcdf-4)

  # These options were removed in version 4.1.3 _sigh_
  if ( ${NetCDF_VERSION} VERSION_LESS 4.1.3 )
    append_set(NetCDF_netcdf4_opts
                        --with-hdf5=${HDF5_install_dir}
                        --with-zlib=${ZLIB_install_dir})
  endif()  

else()   
  set(NetCDF_netcdf4_opts --disable-netcdf-4)
endif() 

# Build CPPFLAGS string. Pick up the CMAKE_BUILD_TYPE flags
include(BuildWhitespaceString)
build_whitespace_string(netcdf_cppflags 
                        -I${TPL_INSTALL_PREFIX}/include ${Amanzi_COMMON_CFLAGS} )

# --- Add external project build and tie to the ZLIB build target
ExternalProject_Add(${NetCDF_BUILD_TARGET}
                    DEPENDS   ${NetCDF_PACKAGE_DEPENDS}             # Package dependency target
                    TMP_DIR   ${NetCDF_tmp_dir}                     # Temporary files directory
                    STAMP_DIR ${NetCDF_stamp_dir}                   # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}              # Download directory
                    URL          ${NetCDF_URL}                    # URL may be a web site OR a local file
                    URL_MD5      ${NetCDF_MD5_SUM}                # md5sum of the archive file
		    # -- Patch 
                    PATCH_COMMAND ${NetCDF_PATCH_COMMAND}
                    # -- Configure
                    SOURCE_DIR       ${NetCDF_source_dir}           # Source directory
                    CONFIGURE_COMMAND
                                    <SOURCE_DIR>/configure
                                                --prefix=<INSTALL_DIR>
                                                --disable-examples
                                                ${NetCDF_netcdf4_opts} 
                                                --disable-dap
                                                --disable-shared
                                                CC=${CMAKE_C_COMPILER}
                                                CXX=${CMAKE_CXX_COMPILER}
                                                FC=${CMAKE_Fortran_COMPILER}
                                                CPPFLAGS=${netcdf_cppflags}
                                                LDFLAGS=-L<INSTALL_DIR>/lib
                    # -- Build
                    BINARY_DIR        ${NetCDF_build_dir}           # Build directory 
                    BUILD_COMMAND     $(MAKE)                     # $(MAKE) enables parallel builds through make
                    BUILD_IN_SOURCE   ${NetCDF_BUILD_IN_SOURCE}     # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}        # Install directory
                    # -- Output control
                    ${NetCDF_logging_args})
