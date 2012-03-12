#  -*- mode: cmake -*-

#
# Build TPL:  HDF5 
#    
define_external_project(HDF5 TARGET hdf5 DEPENDS ZLIB )

# Use the common cflags
include(BuildWhitespaceString)
build_whitespace_string(hdf5_cppflags
                       -I${TPL_INSTALL_PREFIX}/include
                       ${Amanzi_COMMON_CFLAGS})

# ############################################################################ #
# Add external project
# ############################################################################ #
ExternalProject_Add(${HDF5_target}
    DEPENDS zlib
    ${HDF5_ep_directory_args}
    ${HDF5_url_args}
    CONFIGURE_COMMAND
                   <SOURCE_DIR>/configure 
                               --prefix=<INSTALL_DIR>
                               --disable-fortran
                               --disable-cxx
                               --enable-production
                               --enable-largefile
                               --enable-parallel
                               --with-zlib=${ZLIB_DIR}
                               CC=${CMAKE_C_COMPILER}
                               CXX=${CMAKE_CXX_COMPILER}
                               CPPFLAGS=${hdf5_cppflags}
                               LDFLAGS=-L${TPL_INSTALL_PREFIX}/lib
    BUILD_COMMAND $(MAKE)                           
    ${HDF5_logging_opts}                  
)
