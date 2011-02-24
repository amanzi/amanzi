# -*- mode: cmake -*-
# 
# Amanzi Third Party Library (TPL) Definitions
#

# Standard CMake modules see CMAKE_ROOT/Modules
include(CMakeDependentOption)
include(FeatureSummary)

# Amanzi CMake modules see <root source>/tools/cmake
include(PrintVariable)


##############################################################################
# Boost
##############################################################################
find_package( Boost COMPONENTS system filesystem program_options REQUIRED)
set_package_info(Boost
                 "C++ Extension library"
                 "http://www.boost.org"
                 "Required by the MPC")


##############################################################################
# HDF5 - http://www.hdfgroup.org/HDF5/
##############################################################################
#find_package(HDF5 REQUIRED)
#if ( NOT HDF5_IS_PARALLEL ) 
#    message(FATAL_ERROR "The HDF5 installation found in ${HDF5_DIR} is not\n"
#                        "a parallel build. Please re-run cmake and define\n"
#                        "a HDF5 installation that is parallel.\n")
#endif(NOT HDF5_IS_PARALLEL)
set_package_info(HDF5
                "I/O library that creates HDF5 formatted files"
                "http://www.hdfgroup.org/HDF5"
                "Required library for several components in Amanzi"
                )


##############################################################################
# NetCDF - http://www.unidata.ucar.edu/software/netcdf/
##############################################################################
#find_package(NetCDF REQUIRED)
set_package_info(NetCDF
                 "Network Common Data Format (NetCDF)"
                 "http://www.unidata.ucar.edu/software/netcdf/"
                 "Required by ExodusII library")

##############################################################################
# Exodus II -http://sourceforge.net/projects/exodusii
##############################################################################
find_package(ExodusII REQUIRED)
set_package_info(ExodusII
                 "File format library. Originated from Sandia."
                 "http://sourceforge.net/projects/exodusii/"
                 "Required by all the mesh frameworks to read mesh files")


# Enabled TPLs

# Mesh TPLs

# Enable ALL possible mesh frameworks
cmake_dependent_option(ENABLE_ALL_Mesh "Build all Amanzi mesh frameworks" OFF
                       "ENABLE_STK_Mesh;ENABLE_MOAB_Mesh;ENABLE_MSTK" ON) 
set_feature_info(ENABLE_ALL_Mesh
                 "Build all available mesh frameworks"
                  )    

##############################################################################
# STK - Sierra Mesh Tool Kit part of Trilinos
##############################################################################
option(ENABLE_STK_Mesh  "Build Amanzi with the STK mesh framework" ON)
set_feature_info(ENABLE_STK_Mesh
                 "Sierra Mesh Tool Kit (STK Mesh) a Trilinos package"
                 "http://trilinos.sandia.gov/packages/stk"
                 )

##############################################################################
# MOAB - svn co https://svn.mcs.anl.gov/repos/ITAPS/MOAB/trunk MOAB
##############################################################################
option(ENABLE_MOAB_Mesh "Build Amanzi with the MOAB mesh framework" OFF)
set_package_info(ENABLE_MOAB_Mesh
                 "A Mesh-Oriented datABase"
                 "http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB"
                 )

##############################################################################
# MSTK - https://software.lanl.gov/MeshTools/trac/raw-attachment/wiki/WikiStart/mstk-1.80.tar.gz
##############################################################################
cmake_dependent_option(ENABLE_MSTK_Mesh "Build Amanzi with the MSTK mesh framework" OFF
                       "ENABLE_METIS" ON )
set_feature_info(ENABLE_MSTK_Mesh
                 "A mesh framework"
                 "https://software.lanl.gov/MeshTools/trac/wiki"
                 "When enabled will also enable METIS"
                 )

##############################################################################
# METIS - http://glaros.dtc.umn.edu/gkhome/metis/metis/download
##############################################################################
option(ENABLE_METIS "Mesh partitioning library" OFF)
set_feature_info(ENABLE_METIS
                 "Mesh partitioning library"
                 "http://glaros.dtc.umn.edu/gkhome/metis/metis/overview"
                 "Required by MSTK"
                 )




##############################################################################
# CGNS - http://www.cgns.sourceforge.net/
##############################################################################
option(ENABLE_CGNS  "Build Amanzi output library with CGNS" OFF)
set_feature_info(ENABLE_CGNS
                 "CFD General Notation System"
                 "http://cgns.sourceforge.net"
                 "Required to produce VisIt files"
                 )
#if (ENABLE_CGNS)
#    find_package(CGNS REQUIRED)
#endif() 

##############################################################################
# UnitTest++ - http://unittest-cpp.sourceforge.net/
##############################################################################
option(ENABLE_UnitTest "Build Amanzi unit tests. Requires UnitTest++" ON)
set_feature_info(ENABLE_UnitTest
                 "C++ unit test framework"
                 "http://unittest-cpp.sourceforge.net"
                 "Required to run test suite"
                 )








