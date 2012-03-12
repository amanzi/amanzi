#  -*- mode: cmake -*-

#
# Build TPL:  HDF5 
#    
define_external_project(HDF5 TARGET hdf5 DEPENDS ZLIB )

# ############################################################################ #
# Create the patch shell script file
# ############################################################################ #
set(HDF5_sh_patch "${HDF5_prefix_dir}/hdf5-patch-step.sh")
configure_file(${SuperBuild_BUILD_FILES_DIR}/hdf5-patch-step.sh.in
               ${HDF5_sh_patch}
               @ONLY)

# ############################################################################ #
# Add external project
# ############################################################################ #
ExternalProject_Add(${HDF5_target}
    DEPENDS zlib
    ${HDF5_ep_directory_args}
    ${HDF5_url_args}
    PATCH_COMMAND sh ${HDF5_sh_patch}
    CMAKE_ARGS
             -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER}
             -D CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
             -D CMAKE_INSTALL_PREFIX=<INSTALL_DIR>
             -D HDF5_ENABLE_PARALLEL:Bool=TRUE
             -D HDF5_ENABLE_Z_LIB_SUPPORT:bool=TRUE
             -D ZLIB_INCLUDE_DIR=${ZLIB_install_dir}/include
             -D ZLIB_LIBRARY=${ZLIB_install_dir}/lib/libz.a
             -D HDF5_BUILD_HL_LIB:bool=TRUE
             -D HDF5_BUILD_TOOLS:bool=TRUE
    ${HDF5_logging_opts}                  
)

