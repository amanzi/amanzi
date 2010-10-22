SET(CMAKE_CXX_COMPILER mpicxx)
SET(CMAKE_C_COMPILER mpicc)

project (AMANZI_Mesh)
enable_testing()


cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE debug)

set(AMANZI_Mesh_VERSION_MAJOR 0)
set(AMANZI_Mesh_VERSION_MINOR 1)


set(HOME "/home/rao")


##############################################################################
# Curl library - http://www.libcurl.org
##############################################################################

set(Curl_DIR ${HOME}/packages/cmake-2.8.2/Utilities/cmcurl)
find_library(CURL_LIBRARY curl 
             /usr/lib
             /usr/local/lib
             ${Curl_DIR}/lib 
             ${Curl_LIB_DIR})
if (NOT CURL_LIBRARY)
  message(FATAL_ERROR "Error: could not find Curl library"
          "Try defining Curl_DIR: (-D Curl_DIR:FILEPATH=-<Curl_install_prefix>"
          "\nor define Curl_LIB_DIR"
         )
endif()

set(CURL_INCLUDE_DIR ${Curl_DIR}/curl)
find_path(CURL_INCLUDES curl.h 
          /usr/include/curl 
          /usr/local/include/curl 
          ${Curl_DIR}/include/curl 
          ${Curl_DIR}/include 
          ${CURL_INCLUDE_DIR})
if(NOT CURL_INCLUDES)
  message(FATAL_ERROR "Error: could not find Curl includes\n"
          "Try defining Curl_DIR: (-D Curl_DIR:FILEPATH=-<Curl_install_prefix>"
          "\nor define Curl_INCLUDE_DIR"
         )
endif()

add_library(curl SHARED IMPORTED)
set_property(TARGET curl PROPERTY IMPORTED_LOCATION ${CURL_LIBRARY})




##############################################################################
# MPI Message Passing Interface 
# OpenMPI - http:://www.open-mpi.org
# MPICH - http://www.mcs.anl.gov/research/projects/mpich2/
##############################################################################

set(MPI_DIR "/usr/lib/openmpi")
find_library(MPI_LIBRARY mpi
             /usr/lib
             /usr/local/lib
             ${MPI_DIR}/lib
             ${MPI_LIB_DIR})
if(NOT MPI_LIBRARY)
  message(FATAL_ERROR "Error: could not find MPI library\n"
          "Try defining MPI_DIR: (-D MPI_DIR:FILEPATH=<MPI_install_prefix>)"
          "\nor defining MPI_LIB_DIR")
endif()

find_path(MPI_INCLUDES mpi.h 
          /usr/include
          /usr/local/include
          ${MPI_DIR}/include
          ${MPI_INCLUDE_DIR})
if(NOT MPI_INCLUDES)
  message(FATAL_ERROR "Error: could not find MPI includes\n"
          "Try defining MPI_DIR: (-D MPI_DIR:FILEPATH=<MPI_install_prefix>)"
          "\nor defining MPI_INCLUDE_DIR")
endif()

add_library(mpi STATIC IMPORTED)
set_property(TARGET mpi PROPERTY IMPORTED_LOCATION ${MPI_LIBRARY})
set(MPI_INCLUDE_DIR ${MPI_DIR}/include)



##############################################################################
# NetCDF - http://www.unidata.ucar.edu/software/netcdf/
##############################################################################

set (NetCDF_DIR /usr)
find_library(NETCDF_LIBRARY netcdf
             ${NetCDF_DIR}/lib
	     ${NetCDF_LIB_DIR}
             /usr/lib
	     /usr/local/lib
	     NO_DEFAULT_PATH
	     )
 if(NOT NETCDF_LIBRARY)
  message(FATAL_ERROR "Error: could not find netcdf library\n"
          "Try defining NetCDF_DIR: (-D NetCDF_DIR:FILEPATH=<NetCDF_install_prefix>)"
          "\nor defining NetCDF_LIB_DIR")
endif()

find_library(NETCDF_C++_LIBRARY netcdf_c++
             ${NetCDF_DIR}/lib
	     ${NetCDF_LIB_DIR}
             /usr/lib
	     /usr/local/lib
	     NO_DEFAULT_PATH)
 if(NOT NETCDF_C++_LIBRARY)
  message(FATAL_ERROR "Error: could not find netcdf_c++ library\n"
          "Try defining NetCDF_DIR: (-D NetCDF_DIR:FILEPATH=<NetCDF_install_prefix>)"
          "\nor defining NetCDF_LIB_DIR")
endif()

find_path(NETCDF_INCLUDES netcdf.h 
          ${NetCDF_DIR}/include
          ${NetCDF_INCLUDE_DIR}
          /usr/include
          /usr/local/include
	  NO_DEFAULT_PATH)
if(NOT NETCDF_INCLUDES)
  message(FATAL_ERROR "Error: could not find NetCDF includes\n"
          "Try defining NetCDF_DIR: (-D NetCDF_DIR:FILEPATH=<NetCDF_install_prefix>)"
          "\nor defining NetCDF_INCLUDE_DIR")
endif()

#
# DID NOT SEEM TO NEED THIS - but left uncommented in checked in version
#
add_library(netcdf SHARED IMPORTED)
set_property(TARGET netcdf PROPERTY IMPORTED_LOCATION ${NETCDF_LIBRARY})

add_library(netcdf_c++ SHARED IMPORTED)
set_property(TARGET netcdf_c++ PROPERTY IMPORTED_LOCATION ${NETCDF_C++_LIBRARY})

#target_link_libraries(netcdf curl)   # See note under ExodusII.



##############################################################################
# HDF5 - http://www.hdfgroup.org/HDF5/
##############################################################################

find_library(HDF5_LIBRARY hdf5
             /usr/lib
             /usr/local/lib
             ${HDF5_DIR}/lib
             ${HDF5_LIB_DIR})
if(NOT HDF5_LIBRARY)
  message(FATAL_ERROR "Error: could not find HDF5 library\n"
          "Try defining HDF5_DIR: (-D HDF5_DIR:FILEPATH=<HDF5_install_prefix>)"
          "\nor defining HDF5_LIB_DIR")
endif()

find_path(HDF5_INCLUDES hdf5.h 
          /usr/include
          /usr/local/include
          ${HDF5_DIR}/include
          ${HDF5_INCLUDE_DIR})
if(NOT HDF5_INCLUDES)
  message(FATAL_ERROR "Error: could not find HDF5 includes\n"
          "Try defining HDF5_DIR: (-D HDF5_DIR:FILEPATH=<HDF5_install_prefix>)"
          "\nor defining HDF5_INCLUDE_DIR")
endif()

add_library(hdf5 STATIC IMPORTED)
set_property(TARGET hdf5 PROPERTY IMPORTED_LOCATION ${HDF5_LIBRARY})




##############################################################################
# ExodusII - http://sourceforge.net/projects/exodusii/ or 
# get it as part of Trilinos - http:://trilinos.sandia.gov
##############################################################################

set(ExodusII_DIR "${HOME}/packages/ACCESS/libraries/exodusii")
find_library(EXODUSII_LIBRARY exoIIv2c
             ${ExodusII_DIR}/cbind/src
	     ${ExodusII_DIR}/lib
	     ${EXODUS_LIB_DIR})
if(NOT EXODUSII_LIBRARY)
  message(FATAL_ERROR "Error: could not find ExodusII library\n"
          "Try defining ExodusII_DIR: (-D ExodusII_DIR:FILEPATH=<ExodusII_install_prefix>)\n"
          "Or ExodusII_LIB_DIR")
endif()

find_path(EXODUSII_INCLUDES exodusII.h 
          ${ExodusII_DIR}/cbind/include
	  ${ExodusII_DIR}/include
          ${ExodusII_INCLUDE_DIR})
if(NOT EXODUSII_INCLUDES)
  message(FATAL_ERROR "Error: could not find ExodusII includes\n"
          "Try defining ExodusII_DIR: (-D ExodusII_DIR:FILEPATH=<ExodusII_install_prefix>)\n"
          "Or ExodusII_INCLUDE_DIR")
endif()

add_library(exodusii STATIC IMPORTED)
set_property(TARGET exodusii PROPERTY IMPORTED_LOCATION ${EXODUSII_LIBRARY})

# Can't add link dependencies to imported libraries. See http://www.vtk.org/Bug/view.php?id=10395
# target_link_libraries(exodusii netcdf)

# Instead, we create a list of all the libraries required for ExodusII to link.
set(EXODUS_LIBS_ALL exodusii netcdf curl)




##############################################################################
# MOAB - svn co https://svn.mcs.anl.gov/repos/ITAPS/MOAB/trunk MOAB
##############################################################################

set (MOAB_DIR "${HOME}/packages/MOAB")
set (MOAB_LIB_DIR ${MOAB_DIR}/lib)
find_library(MOAB_LIBRARY MOAB
             /usr/lib
             /usr/local/lib
             ${MOAB_LIB_DIR})
if(NOT MOAB_LIBRARY)
  message(FATAL_ERROR "Error: could not find MOAB library\n"
          "Try defining MOAB_DIR: (-D MOAB_DIR:FILEPATH=<MOAB_install_prefix>)"
          "\nor defining MOAB_LIB_DIR")
endif()

set (MOAB_INCLUDE_DIR ${MOAB_DIR}/include)
find_path(MOAB_INCLUDES MBCore.hpp 
          /usr/include
          /usr/local/include
          ${MOAB_INCLUDE_DIR})
if(NOT MOAB_INCLUDES)
  message(FATAL_ERROR "Error: could not find MOAB includes\n"
          "Try defining MOAB_DIR: (-D MOAB_DIR:FILEPATH=<MOAB_install_prefix>)"
          "\nor defining MOAB_INCLUDE_DIR")
endif()

add_library(moab STATIC IMPORTED)
set_property(TARGET moab PROPERTY IMPORTED_LOCATION ${MOAB_LIBRARY})

# create a list of all the libraries required for MOAB to link.
set(MOAB_LIBS_ALL MOAB netcdf netcdf_c++ hdf5)




##############################################################################
# UnitTest++ - http://unittest-cpp.sourceforge.net/
##############################################################################

set(UnitTest_DIR "${HOME}/packages/UnitTest++")
if (NOT DEFINED UnitTest_DIR)
  message(WARNING "No UnitTest++ directory defined. To compile the unit tests add:
 -D UnitTest_DIR::FILEPATH=<UnitTest++_install_prefix>")
else()
  find_library(UNITTEST_LIBRARY unittest++ "${UnitTest_DIR}")
  if (NOT UNITTEST_LIBRARY)
    message(FATAL_ERROR "Error: could not find unittest++ library")
  endif()
  
  find_path(UNITTEST_INCLUDES "UnitTest++.h" "${UnitTest_DIR}/src" "${UnitTest_DIR}/include/unittest++")
  if (NOT UNITTEST_INCLUDES)
    message(FATAL_ERROR "Error: could not find UnitTest++ includes")
  else()
    include_directories(${UNITTEST_INCLUDES})
  endif()

  set(BUILD_TESTS 1)
  add_library(unittest++ STATIC IMPORTED)
  set_property(TARGET unittest++ PROPERTY IMPORTED_LOCATION ${UNITTEST_LIBRARY})

endif()


##############################################################################
#Trilinos
##############################################################################
set (Trilinos_DIR "${HOME}/packages/trilinos/trilinos-10.4.0")
if(NOT DEFINED Trilinos_DIR)
endif()

find_package(Trilinos REQUIRED PATHS ${Trilinos_DIR}/include)
if(NOT Trilinos_FOUND)
  message(FATAL_ERROR "Error: could not find Trilinos\n"
          "Try defining Trilinos_DIR: (-D Trilinos_DIR:FILEPATH=<Trilinos_install_prefix>)")
endif()
list(APPEND Trilinos_LIBRARIES m lapack blas)
list(APPEND Trilinos_DIRS /usr/lib)
list(APPEND Trilinos_INCLUDE_DIRS /usr/include)




include_directories(${PROJECT_BINARY_DIR})

add_subdirectory(dbc)
add_subdirectory(exodus)
add_subdirectory(mesh_data)
add_subdirectory(stk_mesh)
add_subdirectory(simple_mesh)
add_subdirectory(moab_mesh)
add_subdirectory(mpc)
add_subdirectory(chemistry)
add_subdirectory(transport)
add_subdirectory(boundary_conditions)
add_subdirectory(flow)
add_subdirectory(gmv)
add_subdirectory(bin)
