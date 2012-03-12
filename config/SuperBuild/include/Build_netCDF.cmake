#  -*- mode: cmake -*-

#
# Build TPL: NETCDF 
#    
define_external_project(NETCDF 
                        TARGET netcdf
                        DEPENDS ${CURL_target} ${HDF5_target}
            )

# ########################################################################### #
# Build the patch command
# ########################################################################### #

# Need Perl to patch the files
find_package(Perl)
if (NOT PERL_FOUND)
  message(FATAL_ERROR "Can not locate Perl. Unable to patch and build netCDF")
endif()


# Configure the bash patch script
set(NETCDF_sh_patch ${NETCDF_prefix_dir}/netcdf-patch-step.sh)
configure_file(${SuperBuild_BUILD_FILES_DIR}/netcdf-patch-step.sh.in
               ${NETCDF_sh_patch}
       @ONLY)

# Configure the CMake command file
set(NETCDF_cmake_patch ${NETCDF_prefix_dir}/netcdf-patch-step.cmake)
configure_file(${SuperBuild_BUILD_FILES_DIR}/netcdf-patch-step.cmake.in
               ${NETCDF_cmake_patch}
       @ONLY)
set(NETCDF_PATCH_COMMAND ${CMAKE_COMMAND} -P ${NETCDF_cmake_patch})     

# ########################################################################### #
# Build the configure command
# ########################################################################### #

option(ENABLE_NETCDF4 "Enable netCDF4 build" TRUE)
message(STATUS "NETCDF4 enabled == ${ENABLE_NETCDF4}")
set(NETCDF_netcdf4_opts)
if (ENABLE_NETCDF4)

  append_set(NETCDF_netcdf4_opts --enable-netcdf-4)

  # These options were removed in version 4.1.3
  if ( ${NETCDF_VERSION} VERSION_LESS 4.1.3 )
    append_set(NETCDF_netcdf4_opts
                        --with-hdf5=${HDF5_install_dir}
                        --with-zlib=${ZLIB_install_dir})
  endif()  

else()   
  set(NETCDF_netcdf4_opts --disable-netcdf-4)
endif() 

# Build CPPFLAGS string. Pick up the CMAKE_BUILD_TYPE flags
include(BuildWhitespaceString)
build_whitespace_string(netcdf_cppflags 
                        -I${TPL_INSTALL_PREFIX}/include ${Amanzi_COMMON_CFLAGS} )

# Add the external project

# Define the compilers through environment variables
# configure script does not pick up these definitions
# from the configure command line
ExternalProject_Add(${NETCDF_target}
    DEPENDS curl hdf5
    ${NETCDF_ep_directory_args}
    ${NETCDF_url_args}
    PATCH_COMMAND ${NETCDF_PATCH_COMMAND}
    CONFIGURE_COMMAND <SOURCE_DIR>/configure
                                   --prefix=<INSTALL_DIR>
                                   --disable-examples
                                   ${NETCDF_netcdf4_opts} 
                                   --disable-dap
                                   --disable-shared
                                   CC=${CMAKE_C_COMPILER}
                                   CXX=${CMAKE_CXX_COMPILER}
                                   FC=${CMAKE_Fortran_COMPILER}
                                   CPPFLAGS=${netcdf_cppflags}
                                   LDFLAGS=-L<INSTALL_DIR>/lib
   BUILD_COMMAND $(MAKE)   				   
   ${NETCDF_logging_opts}   
)
