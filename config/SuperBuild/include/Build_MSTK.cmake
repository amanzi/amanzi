#  -*- mode: cmake -*-

#
# Build TPL: MSTK 
#    
# --- Define all the directories and common external project flags
define_external_project_args(MSTK
                             TARGET mstk
                             DEPENDS HDF5 NetCDF ExodusII METIS)


# --- Define the configure parameters

# Compile flags
set(mstk_cflags_list -I${TPL_INSTALL_PREFIX}/include ${Amanzi_COMMON_CFLAGS})
build_whitespace_string(mstk_cflags ${mstk_cflags_list})

# The CMake cache args
set(MSTK_CMAKE_CACHE_ARGS
                    ${Amanzi_CMAKE_C_COMPILER_ARGS}
                    -DCMAKE_C_FLAGS:STRING=${mstk_cflags}
                    -DCMAKE_EXE_LINKER_FLAGS:STRING=-L${TPL_INSTALL_PREFIX}/lib
                    -DENABLE_PARALLEL:BOOL=TRUE
                    -DENABLE_ExodusII:BOOL=TRUE
                    -DENABLE_ZOLTAN:BOOL=TRUE
                    -DHDF5_DIR:PATH=${TPL_INSTALL_PREFIX}
                    -DNetCDF_DIR:PATH=${TPL_INSTALL_PREFIX} 
                    -DExodusII_DIR:PATH=${TPL_INSTALL_PREFIX} 
                    -DZOLTAN_DIR:PATH=${TPL_INSTALL_PREFIX}
                    -DMetis_DIR:PATH=${TPL_INSTALL_PREFIX} 
                    -DMETIS_LIB_DIR:PATH=${TPL_INSTALL_PREFIX}/lib 
                    -DMETIS_LIBRARY:PATH=${METIS_LIBRARY}
                    -DMetis_INCLUDE_DIR:PATH=${TPL_INSTALL_PREFIX}/include 
                    -DMETIS_INCLUDE_DIRS:PATH=${TPL_INSTALL_PREFIX}/include
                    -DENABLE_Tests:BOOL=FALSE
                    -DINSTALL_DIR:PATH=<INSTALL_DIR>
                    -DINSTALL_ADD_VERSION:BOOL=FALSE)

# --- Add external project build and tie to the MSTK build target
ExternalProject_Add(${MSTK_BUILD_TARGET}
                    DEPENDS   ${MSTK_PACKAGE_DEPENDS}             # Package dependency target
                    TMP_DIR   ${MSTK_tmp_dir}                     # Temporary files directory
                    STAMP_DIR ${MSTK_stamp_dir}                   # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}              # Download directory
                    URL          ${MSTK_URL}                      # URL may be a web site OR a local file
                    URL_MD5      ${MSTK_MD5_SUM}                  # md5sum of the archive file
                    # -- Configure
                    SOURCE_DIR       ${MSTK_source_dir}           # Source directory
                    CMAKE_CACHE_ARGS ${MSTK_CMAKE_CACHE_ARGS}
                    # -- Build
                    BINARY_DIR        ${MSTK_build_dir}           # Build directory 
                    BUILD_COMMAND     $(MAKE)                     # $(MAKE) enables parallel builds through make
                    BUILD_IN_SOURCE   ${MSTK_BUILD_IN_SOURCE}     # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}        # Install directory
                    # -- Output control
                    ${MSTK_logging_args})


# MSTK include and library install path
global_set(MSTK_INCLUDE_DIR "${TPL_INSTALL_PREFIX}/include")
global_set(MSTK_ARCHOS "${CMAKE_SYSTEM_PROCESSOR}_${CMAKE_SYSTEM_NAME}")
global_set(MSTK_LIBRARY_DIR "${TPL_INSTALL_PREFIX}/lib/${MSTK_ARCHOS}")
