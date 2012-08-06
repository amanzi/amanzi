#  -*- mode: cmake -*-

#
# Build TPL: NetCDF 
# 

# --- Define all the directories and common external project flags
define_external_project_args(NetCDF 
                             TARGET netcdf
                             DEPENDS ${MPI_PROJECT} CURL
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
                        --with-netcdf=${HDF5_install_dir}
                        --with-zlib=${ZLIB_install_dir})
  endif()  

else()   
  set(NetCDF_netcdf4_opts --disable-netcdf-4)
endif() 

# Build compiler flag strings for C, C++ and Fortran
include(BuildWhitespaceString)
build_whitespace_string(netcdf_cflags 
                        -I${TPL_INSTALL_PREFIX}/include ${Amanzi_COMMON_CFLAGS} )

build_whitespace_string(netcdf_cxxflags 
                        -I${TPL_INSTALL_PREFIX}/include ${Amanzi_COMMON_CXXFLAGS} )

set(cpp_flags_list
    -I${TPL_INSTALL_PREFIX}/include
    ${Amanzi_COMMON_CFLAGS}
    ${Amanzi_COMMON_CXXFLAGS})
list(REMOVE_DUPLICATES cpp_flags_list)
build_whitespace_string(netcdf_cppflags ${cpp_flags_list})

build_whitespace_string(netcdf_fcflags 
                        ${Amanzi_COMMON_FCFLAGS} )

# Add MPI C libraries 
if ( ( NOT BUILD_MPI) AND ( NOT MPI_WRAPPERS_IN_USE ) AND (MPI_C_LIBRARIES) )
  build_whitespace_string(netcdf_ldflags -L${TPL_INSTALL_PREFIX}/lib ${MPI_C_LIBRARIES} ${CMAKE_EXE_LINKER_FLAGS})
else()
  build_whitespace_string(netcdf_ldflags -L${TPL_INSTALL_PREFIX}/lib ${CMAKE_EXE_LINKER_FLAGS})
endif()  


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
                                                --disable-fortran
                                                --disable-f90
                                                --disable-f77
                                                --disable-fortran-compiler-check
                                                CC=${CMAKE_C_COMPILER_USE}
                                                CFLAGS=${netcdf_cflags}
                                                CXX=${CMAKE_CXX_COMPILER_USE}
                                                CXXFLAGS=${netcdf_cxxflags}
                                                CPPFLAGS=${netcdf_cppflags}
                                                LDFLAGS=${netcdf_ldflags}
                    # -- Build
                    BINARY_DIR        ${NetCDF_build_dir}           # Build directory 
                    BUILD_COMMAND     $(MAKE)                     # $(MAKE) enables parallel builds through make
                    BUILD_IN_SOURCE   ${NetCDF_BUILD_IN_SOURCE}     # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}        # Install directory
                    # -- Output control
                    ${NetCDF_logging_args})
