#  -*- mode: cmake -*-

#
# Build TPL: SEACAS 
#    
# --- Define all the directories and common external project flags

# SEACAS does not call MPI directly, however HDF5 requires
# MPI and to resolve links we need MPI compile wrappers.
define_external_project_args(SEACAS
                             TARGET seacas
                             DEPENDS HDF5 NetCDF)


# --- Define the configure parameters

# Compile flags
set(seacas_cflags_list -I${TPL_INSTALL_PREFIX}/include ${Amanzi_COMMON_CFLAGS})
build_whitespace_string(seacas_cflags ${seacas_cflags_list})

set(seacas_cxxflags_list -I${TPL_INSTALL_PREFIX}/include ${Amanzi_COMMON_CXXFLAGS})
build_whitespace_string(seacas_cflags ${seacas_cxxflags_list})

set(seacas_fcflags_list -I${TPL_INSTALL_PREFIX}/include ${Amanzi_COMMON_FCFLAGS})
build_whitespace_string(seacas_fcflags ${seacas_fcflags_list})

# Build the NetCDF libraries string
include(BuildLibraryName)
build_library_name(netcdf seacas_netcdf_library STATIC APPEND_PATH ${TPL_INSTALL_PREFIX}/lib)
build_library_name(hdf5_hl seacas_hdf5_hl_library STATIC APPEND_PATH ${TPL_INSTALL_PREFIX}/lib)
build_library_name(hdf5 seacas_hdf5_library STATIC APPEND_PATH ${TPL_INSTALL_PREFIX}/lib)
build_library_name(z seacas_z_library STATIC APPEND_PATH ${TPL_INSTALL_PREFIX}/lib)
set(seacas_netcdf_libraries
       ${seacas_netcdf_library}
       ${seacas_hdf5_hl_library}
       ${seacas_hdf5_library}
       ${seacas_z_library})

# The CMake cache args
set(SEACAS_CMAKE_CACHE_ARGS
                    -DCMAKE_INSTALL_PREFIX:FILEPATH=<INSTALL_DIR>
                    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
                    ${Amanzi_CMAKE_C_COMPILER_ARGS}
                    -DCMAKE_C_COMPILER:FILEPATH=${MPI_C_COMPILER}
                    ${Amanzi_CMAKE_CXX_COMPILER_ARGS}
                    -DCMAKE_CXX_COMPILER:FILEPATH=${MPI_CXX_COMPILER}
                    ${Amanzi_CMAKE_Fortran_COMPILER_ARGS}
                    -DCMAKE_Fortran_COMPILER:FILEPATH=${MPI_Fortran_COMPILER}
                    -DCMAKE_EXE_LINKER_FLAGS:STRING=-L${TPL_INSTALL_PREFIX}/lib
                    -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=FALSE
                    -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=FALSE
                    -DTrilinos_ENABLE_SEACAS:BOOL=TRUE
                    -DTPL_Netcdf_LIBRARIES:STRING=${seacas_netcdf_libraries}
                    -DNetcdf_INCLUDE_DIRS:STRING=${TPL_INSTALL_PREFIX}/include
                    )

# --- Add external project build and tie to the SEACAS build target
ExternalProject_Add(${SEACAS_BUILD_TARGET}
                    DEPENDS   ${SEACAS_PACKAGE_DEPENDS}             # Package dependency target
                    TMP_DIR   ${SEACAS_tmp_dir}                     # Temporary files directory
                    STAMP_DIR ${SEACAS_stamp_dir}                   # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}                # Download directory
                    URL          ${SEACAS_URL}                      # URL may be a web site OR a local file
                    URL_MD5      ${SEACAS_MD5_SUM}                  # md5sum of the archive file
                    # -- Configure
                    SOURCE_DIR       ${SEACAS_source_dir}           # Source directory
                    CMAKE_CACHE_ARGS ${SEACAS_CMAKE_CACHE_ARGS}
                    # -- Build
                    BINARY_DIR        ${SEACAS_build_dir}           # Build directory 
                    BUILD_COMMAND     $(MAKE)                       # $(MAKE) enables parallel builds through make
                    BUILD_IN_SOURCE   ${SEACAS_BUILD_IN_SOURCE}     # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}/SEACAS   # Install directory, NOT in the usual place!
                    # -- Output control
                    ${SEACAS_logging_args})
