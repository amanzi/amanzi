#  -*- mode: cmake -*-

#
# Build TPL: MSTK 
#    
define_external_project(MSTK 
                        TARGET mstk
                        DEPENDS metis hdf5 netcdf exodusii) 


# Define the compile flags
set(mstk_cflags_list -I${TPL_INSTALL_PREFIX}/include ${Amanzi_COMMON_CFLAGS})
build_whitespace_string(mstk_cflags ${mstk_cflags_list})
                            
# ########################################################################### # 
# Add External Project 
# ########################################################################### # 
ExternalProject_Add(${MSTK_target}
    DEPENDS metis hdf5 netcdf exodusii
    ${MSTK_ep_directory_args}
    ${MSTK_url_args}
    CMAKE_ARGS
                    ${Amanzi_CMAKE_COMPILER_ARGS}
                    -DCMAKE_C_FLAGS=${mstk_cflags}
                    -DCMAKE_EXE_LINKER_FLAGS:STRING=-L${TPL_INSTALL_PREFIX}/lib
                    -DENABLE_PARALLEL=yes 
                    -DENABLE_ExodusII=yes 
                    -DHDF5_DIR=${HDF5_install_dir}
                    -DNetCDF_DIR:FILEPATH=${NETCDF_install_dir} 
                    -DExodusII_DIR:FILEPATH=${ExodusII_install_dir} 
                    -DMetis_DIR:FILEPATH=${METIS_install_dir} 
                    -DMETIS_LIB_DIR:FILEPATH=${METIS_install_dir}/lib 
                    -DMETIS_LIBRARY=${METIS_LIBRARY}
                    -DMetis_INCLUDE_DIR:FILEPATH=${METIS_install_dir}/include 
                    -DMETIS_INCLUDES=${TPL_INSTALL_PREFIX}/include
                    -DENABLE_Tests=no 
                    -DINSTALL_DIR:PATH=<INSTALL_DIR>
                    -DINSTALL_ADD_VERSION=no 
    ${MSTK_logging_opts}               
)

# MSTK include and library install path
set(MSTK_INCLUDE_DIR "${MSTK_install_dir}/include")
set(MSTK_ARCHOS "${CMAKE_SYSTEM_PROCESSOR}_${CMAKE_SYSTEM_NAME}")
set(MSTK_LIBRARY_DIR "${MSTK_install_dir}/lib/${MSTK_ARCHOS}")
