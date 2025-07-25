#  -*- mode: cmake -*-

#
# Build TPL: NetCDF 
# 
# --- Define all the directories and common external project flags
define_external_project_args(NetCDF 
                             TARGET netcdf
                             DEPENDS HDF5)

# add version version to the autogenerated tpl_versions.h file
amanzi_tpl_version_write(FILENAME ${TPL_VERSIONS_INCLUDE_FILE}
  PREFIX NetCDF
  VERSION ${NetCDF_VERSION_MAJOR} ${NetCDF_VERSION_MINOR} ${NetCDF_VERSION_PATCH})

# --- Patch the original code
set(NetCDF_patch_file netcdf-cmake.patch
                      netcdf-cmake-rpath.patch)
#                     netcdf-cmake-dl.patch)
#                     netcdf-cmake-namespace.patch)
set(NetCDF_sh_patch ${NetCDF_prefix_dir}/netcdf-patch-step.sh)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/netcdf-patch-step.sh.in
               ${NetCDF_sh_patch}
               @ONLY)

# configure the CMake patch step
set(NetCDF_cmake_patch ${NetCDF_prefix_dir}/netcdf-patch-step.cmake)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/netcdf-patch-step.cmake.in
               ${NetCDF_cmake_patch}
               @ONLY)

# configure the CMake command file
set(NetCDF_PATCH_COMMAND ${CMAKE_COMMAND} -P ${NetCDF_cmake_patch})     

# --- Define the configure command
set(NetCDF_CMAKE_CACHE_ARGS "-DCMAKE_INSTALL_PREFIX:FILEPATH=${TPL_INSTALL_PREFIX}")
list(APPEND NetCDF_CMAKE_CACHE_ARGS "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}")
list(APPEND NetCDF_CMAKE_CACHE_ARGS "-DCMAKE_INSTALL_LIBDIR:FILEPATH=lib")
list(APPEND NetCDF_CMAKE_CACHE_ARGS "-DCMAKE_INSTALL_BINDIR:FILEPATH=bin")
list(APPEND NetCDF_CMAKE_CACHE_ARGS "-DENABLE_DAP:BOOL=FALSE")
list(APPEND NetCDF_CMAKE_CACHE_ARGS "-DENABLE_PARALLEL4:BOOL=TRUE")
list(APPEND NetCDF_CMAKE_CACHE_ARGS "-DHDF5_PARALLEL:BOOL=TRUE")
list(APPEND NetCDF_CMAKE_CACHE_ARGS "-DHDF5_C_LIBRARY:FILEPATH=${HDF5_C_LIBRARY}")
list(APPEND NetCDF_CMAKE_CACHE_ARGS "-DHDF5_HL_LIBRARY:FILEPATH=${HDF5_HL_LIBRARY}")
list(APPEND NetCDF_CMAKE_CACHE_ARGS "-DHDF5_INCLUDE_DIR:PATH=${HDF5_INCLUDE_DIRS}")
list(APPEND NetCDF_CMAKE_CACHE_ARGS "-DHDF5_VERSION:STRING=${HDF5_VERSION}")

# --- Override minimum version
if(CMAKE_MAJOR_VERSION VERSION_EQUAL "4")
  list(APPEND NetCDF_CMAKE_CACHE_ARGS "-DCMAKE_POLICY_VERSION_MINIMUM:STRING=3.5")
endif()

# Default is to build with NetCDF4 which depends on HDF5
option(ENABLE_NetCDF4 "Enable netCDF4 build" TRUE)
if (ENABLE_NetCDF4)
  list(APPEND NetCDF_PACKAGE_DEPENDS ${HDF5_BUILD_TARGET})
  list(APPEND NetCDF_CMAKE_CACHE_ARGS "-DENABLE_NETCDF_4:BOOL=TRUE")
  list(APPEND NetCDF_CMAKE_CACHE_ARGS "-DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARIES}")
endif() 

# share libraries -- disabled by default
list(APPEND NetCDF_CMAKE_CACHE_ARGS "-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}")
if (BUILD_STATIC_LIBS)
  list(APPEND NetCDF_CMAKE_CACHE_ARGS "-DCMAKE_EXE_LINKER_FLAGS:STRING=-ldl")
endif()

# --- Add external project build 
ExternalProject_Add(${NetCDF_BUILD_TARGET}
                    DEPENDS   ${NetCDF_PACKAGE_DEPENDS}   # Package dependency target
                    TMP_DIR   ${NetCDF_tmp_dir}           # Temporary files directory
                    STAMP_DIR ${NetCDF_stamp_dir}         # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}      # Download directory
                    URL          ${NetCDF_URL}            # URL may be a web site OR a local file
                    URL_MD5      ${NetCDF_MD5_SUM}        # md5sum of the archive file
                    DOWNLOAD_NAME ${NetCDF_SAVEAS_FILE}   # file name to store (if not end of URL)
                    # -- Patch 
                    PATCH_COMMAND ${NetCDF_PATCH_COMMAND}
                    # -- Configure
                    SOURCE_DIR       ${NetCDF_source_dir}
                    CMAKE_CACHE_ARGS ${AMANZI_CMAKE_CACHE_ARGS}  # Ensure uniform build
                                     ${NetCDF_CMAKE_CACHE_ARGS}
                                     -DCMAKE_C_FLAGS:STRING=${Amanzi_COMMON_CFLAGS}  # Ensure uniform build
                                     -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                                     -DCMAKE_CXX_FLAGS:STRING=${Amanzi_COMMON_CXXFLAGS}
                                     -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                    # -- Build
                    BINARY_DIR        ${NetCDF_build_dir}       # Build directory 
                    BUILD_COMMAND     $(MAKE)                   # $(MAKE) enables parallel builds through make
                    BUILD_IN_SOURCE   ${NetCDF_BUILD_IN_SOURCE} # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}      
                    # -- Output control
                    ${NetCDF_logging_args})

# --- Useful variables for packages that depend on NetCDF (Trilinos)
include(BuildLibraryName)
build_library_name(netcdf NetCDF_C_LIBRARY APPEND_PATH ${TPL_INSTALL_PREFIX}/lib)

build_library_name(netcdf_c++ NetCDF_CXX_LIBRARY APPEND_PATH ${TPL_INSTALL_PREFIX}/lib)
set(NetCDF_DIR ${TPL_INSTALL_PREFIX})
set(NetCDF_INCLUDE_DIRS ${TPL_INSTALL_PREFIX}/include)
set(NetCDF_C_LIBRARIES ${NetCDF_C_LIBRARY})
if (ENABLE_NetCDF4)
  list(APPEND NetCDF_C_LIBRARIES ${HDF5_LIBRARIES})
  list(APPEND NetCDF_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
  list(REMOVE_DUPLICATES NetCDF_INCLUDE_DIRS)
endif()
  

