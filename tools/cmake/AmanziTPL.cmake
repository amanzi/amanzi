# -*- mode: cmake -*-
# 
# Amanzi Thrid Party Library (TPL) Definitions
#

# Amanzi CMake modules see <root source>/tools/cmake
include(PrintVariable)

# Required TPLs
set(AMANZI_REQ_TPLS "NetCDF" "HDF5" CACHE INTERNAL "List of required TPLs to build Amanzi")

##############################################################################
# HDF5 - http://www.hdfgroup.org/HDF5/
##############################################################################
find_package(HDF5 REQUIRED)
if ( NOT HDF5_IS_PARALLEL ) 
    message(FATAL_ERROR "The HDF5 installation found in ${HDF5_DIR} is not"
                        "a parallel build. Please re-run cmake and define"
                        "a HDF5 installation that is parallel.")
endif(NOT HDF5_IS_PARALLEL)

##############################################################################
# NetCDF - http://www.unidata.ucar.edu/software/netcdf/
##############################################################################
find_package(NetCDF REQUIRED)

##############################################################################
# Exodus II -http://sourceforge.net/projects/exodusii
##############################################################################
find_package(ExodusII REQUIRED)


# Enabled TPLs
option(ENABLE_STK_Mesh  "Build Amanzi with the STK mesh framework" ON)
option(ENABLE_MOAB_Mesh "Build Amanzi with the MOAB mesh framework" OFF)
option(ENABLE_MSTK_Mesh "Build Amanzi with the MSTK mesh framework" OFF)

##############################################################################
# CGNS - http://www.cgns.sourceforge.net/
##############################################################################
option(ENABLE_CGNS      "Build Amanzi output library with CGNS" OFF)
if (ENABLE_CGNS)
    find_package(CGNS REQUIRED)
endif()    




