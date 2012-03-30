#  -*- mode: cmake -*-

#
# Build TPL: MOAB 
# 

# --- Define all the directories and common external project flags
define_external_project_args(MOAB 
                             TARGET moab
                             DEPENDS ZLIB HDF5 NetCDF
                            )

# --- Build common compiler and link flags

# Build compiler flag strings for C
include(BuildWhitespaceString)
build_whitespace_string(moab_cflags 
                        -I${TPL_INSTALL_PREFIX}/include ${Amanzi_COMMON_CFLAGS} )

build_whitespace_string(moab_cxxflags 
                        -I${TPL_INSTALL_PREFIX}/include ${Amanzi_COMMON_CXXFLAGS} )

# Build the LDFLAGS string		      
build_whitespace_string(moab_ldflags
                        -L<INSTALL_DIR>/lib
			-L${TPL_INSTALL_PREFIX}/lib
			-lnetcdf
			-L${TPL_INSTALL_PREFIX}/lib
			-lhdf5_hl
			-lhdf5
			-L${TPL_INSTALL_PREFIX}/lib
			-lz)

# --- Add external project build and tie to the ZLIB build target
ExternalProject_Add(${MOAB_BUILD_TARGET}
                    DEPENDS   ${MOAB_PACKAGE_DEPENDS}             # Package dependency target
                    TMP_DIR   ${MOAB_tmp_dir}                     # Temporary files directory
                    STAMP_DIR ${MOAB_stamp_dir}                   # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}              # Download directory
                    URL          ${MOAB_URL}                    # URL may be a web site OR a local file
                    URL_MD5      ${MOAB_MD5_SUM}                # md5sum of the archive file
                    # -- Configure
                    SOURCE_DIR       ${MOAB_source_dir}           # Source directory
                    CONFIGURE_COMMAND
                                    <SOURCE_DIR>/configure
                                                --prefix=<INSTALL_DIR>
                                                --disable-fortran
						--with-mpi
						--with-hdf5=${TPL_INSTALL_PREFIX}
                                                --with-netcdf=${TPL_INSTALL_PREFIX}
                                                CC=${CMAKE_C_COMPILER}
                                                CFLAGS=${moab_cflags}
                                                CXX=${CMAKE_CXX_COMPILER}
                                                CFLAGS=${moab_cxxflags}
                                                LDFLAGS=${moab_ldflags}
                    # -- Build
                    BINARY_DIR        ${MOAB_build_dir}           # Build directory 
                    BUILD_COMMAND     $(MAKE)                     # $(MAKE) enables parallel builds through make
                    BUILD_IN_SOURCE   ${MOAB_BUILD_IN_SOURCE}     # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}        # Install directory
                    # -- Output control
                    ${MOAB_logging_args})
