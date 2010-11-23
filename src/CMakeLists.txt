# -*- mode: cmake -*-
# Note: This CMakeLists.txt file previously hard-wired the C and C++ compilers 
# via the following commands:
#
#   SET(CMAKE_CXX_COMPILER mpicxx)
#   SET(CMAKE_C_COMPILER mpicc)
#
# The values of these variables should not be set within this file; on some 
# systems no "mpicxx" and "mpicc" wrappers are present.  Instead, cmake 
# should be invoked with arguments such as
#
# cmake \
#   -D CMAKE_CXX_COMPILER:FILEPATH="${ASYNCPE_DIR}/bin/CC" \
#   -D CMAKE_C_COMPILER:FILEPATH="${ASYNCPE_DIR}/bin/cc" \
#
# --RTM
    
project (AMANZI)
enable_testing()


cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE debug)

set(AMANZI_Mesh_VERSION_MAJOR 0)
set(AMANZI_Mesh_VERSION_MINOR 1)

##############################################################################
# Attempting to explicitly find the MPI implementation libraries has been 
# deprecated in favor of using compiler wrappers such as 'mpicc'.  See Trac 
# ticket #240.  --RTM
##############################################################################
# MPI Message Passing Interface 
# OpenMPI - http:://www.open-mpi.org
# MPICH - http://www.mcs.anl.gov/research/projects/mpich2/
##############################################################################
#
#find_library(MPI_LIBRARY
#             NAMES mpi mpich
#             HINTS ${MPI_LIB_DIR}
#                   ${MPI_DIR}/lib
#                   /usr/lib64
#                   /usr/local/lib)
#if(NOT MPI_LIBRARY)
#  message(FATAL_ERROR "Error: could not find MPI library\n"
#          "Try defining MPI_DIR: (-D MPI_DIR:FILEPATH=<MPI_install_prefix>)"
#          "\nor defining MPI_LIB_DIR")
#endif()
#
#find_path(MPI_INCLUDES mpi.h 
#          HINTS ${MPI_INCLUDE_DIR}
#                ${MPI_DIR}/include
#                /usr/include/mpi
#                /usr/include
#                /usr/local/include)
#if(NOT MPI_INCLUDES)
#  message(FATAL_ERROR "Error: could not find MPI includes\n"
#          "Try defining MPI_DIR: (-D MPI_DIR:FILEPATH=<MPI_install_prefix>)"
#          "\nor defining MPI_INCLUDE_DIR")
#endif()
#
#add_library(mpi STATIC IMPORTED)
#set_property(TARGET mpi PROPERTY IMPORTED_LOCATION ${MPI_LIBRARY})
#
#set(MPI_LIBS_ALL mpi mpi_cxx)

##############################################################################
# NetCDF - http://www.unidata.ucar.edu/software/netcdf/
##############################################################################

find_library(NETCDF_LIBRARY netcdf
	     HINTS ${NetCDF_LIB_DIR}
                   ${NetCDF_DIR}/lib
                   /usr/lib
	           /usr/local/lib
		   )
if(NOT NETCDF_LIBRARY)
  message(FATAL_ERROR "Error: could not find netcdf library\n"
          "Try defining NetCDF_DIR: (-D NetCDF_DIR:FILEPATH=<NetCDF_install_prefix>)"
          "\nor defining NetCDF_LIB_DIR")
endif()

find_library(NETCDF_C++_LIBRARY netcdf_c++
             HINTS ${NetCDF_LIB_DIR}
                   ${NetCDF_DIR}/lib
                   /usr/lib
	           /usr/local/lib
		   )
 if(NOT NETCDF_C++_LIBRARY)
  message(FATAL_ERROR "Error: could not find netcdf_c++ library\n"
          "Try defining NetCDF_DIR: (-D NetCDF_DIR:FILEPATH=<NetCDF_install_prefix>)"
          "\nor defining NetCDF_LIB_DIR")
endif()

find_path(NETCDF_INCLUDES netcdf.h 
          HINTS ${NetCDF_INCLUDE_DIR}
                ${NetCDF_DIR}/include
                /usr/include
                /usr/local/include
		)
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
set(NETCDF_LIBS_ALL netcdf curl hdf5 hdf5_hl z mpi mpi_cxx)

##############################################################################
# HDF5 - http://www.hdfgroup.org/HDF5/
##############################################################################

find_library(HDF5_LIBRARY hdf5
             HINTS ${HDF5_LIB_DIR}
                   ${HDF5_DIR}/lib
                   /usr/lib
                   /usr/local/lib
		   )
if(NOT HDF5_LIBRARY)
  message(FATAL_ERROR "Error: could not find HDF5 library\n"
          "Try defining HDF5_DIR: (-D HDF5_DIR:FILEPATH=<HDF5_install_prefix>)"
          "\nor defining HDF5_LIB_DIR")
endif()

find_library(HDF5_HL_LIBRARY hdf5_hl
             HINTS ${HDF5_LIB_DIR}
                   ${HDF5_DIR}/lib
                   /usr/lib
                   /usr/local/lib
		   )
if(NOT HDF5_HL_LIBRARY)
  message(FATAL_ERROR "Error: could not find HDF5_HL library\n"
          "Try defining HDF5_DIR: (-D HDF5_DIR:FILEPATH=<HDF5_install_prefix>)"
          "\nor defining HDF5_LIB_DIR")
endif()

find_path(HDF5_INCLUDES hdf5.h 
          HINTS ${HDF5_INCLUDE_DIR}
                ${HDF5_DIR}/include
                /usr/include
                /usr/local/include
		)
if(NOT HDF5_INCLUDES)
  message(FATAL_ERROR "Error: could not find HDF5 includes\n"
          "Try defining HDF5_DIR: (-D HDF5_DIR:FILEPATH=<HDF5_install_prefix>)"
          "\nor defining HDF5_INCLUDE_DIR")
endif()

add_library(hdf5 STATIC IMPORTED)
set_property(TARGET hdf5 PROPERTY IMPORTED_LOCATION ${HDF5_LIBRARY})

add_library(hdf5_hl STATIC IMPORTED)
set_property(TARGET hdf5_hl PROPERTY IMPORTED_LOCATION ${HDF5_HL_LIBRARY})

set(HDF5_LIBS_ALL hdf5 hdf5_hl z mpi mpi_cxx)



##############################################################################
# ExodusII - http://sourceforge.net/projects/exodusii/ or 
# get it as part of Trilinos - http:://trilinos.sandia.gov
##############################################################################

find_library(EXODUSII_LIBRARY exoIIv2c
	     HINTS ${EXODUS_LIB_DIR}
                   ${ExodusII_DIR}/cbind/src
	           ${ExodusII_DIR}/lib
		   )
if(NOT EXODUSII_LIBRARY)
  message(FATAL_ERROR "Error: could not find ExodusII library\n"
          "Try defining ExodusII_DIR: (-D ExodusII_DIR:FILEPATH=<ExodusII_install_prefix>)\n"
          "Or ExodusII_LIB_DIR")
endif()

find_path(EXODUSII_INCLUDES exodusII.h 
          HINTS ${ExodusII_INCLUDE_DIR}
                ${ExodusII_DIR}/cbind/include
	        ${ExodusII_DIR}/include
		)
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
set(EXODUS_LIBS_ALL exodusii netcdf hdf5 hdf5_hl z curl mpi mpi_cxx) 




##############################################################################
# MOAB - svn co https://svn.mcs.anl.gov/repos/ITAPS/MOAB/trunk MOAB
##############################################################################

find_library(MOAB_LIBRARY MOAB
             HINTS ${MOAB_LIB_DIR}
                   ${MOAB_DIR}/lib
                   /usr/lib
                   /usr/local/lib
		   )
if(NOT MOAB_LIBRARY)
  message(FATAL_ERROR "Error: could not find MOAB library\n"
          "Try defining MOAB_DIR: (-D MOAB_DIR:FILEPATH=<MOAB_install_prefix>)"
          "\nor defining MOAB_LIB_DIR")
endif()

find_path(MOAB_INCLUDES MBCore.hpp 
          HINTS ${MOAB_INCLUDE_DIR}
                ${MOAB_DIR}/include
                /usr/include
                /usr/local/include
		)
if(NOT MOAB_INCLUDES)
  message(FATAL_ERROR "Error: could not find MOAB includes\n"
          "Try defining MOAB_DIR: (-D MOAB_DIR:FILEPATH=<MOAB_install_prefix>)"
          "\nor defining MOAB_INCLUDE_DIR")
endif()

add_library(moab STATIC IMPORTED)
set_property(TARGET moab PROPERTY IMPORTED_LOCATION ${MOAB_LIBRARY})

# create a list of all the libraries required for MOAB to link.
set(MOAB_LIBS_ALL moab netcdf_c++ netcdf hdf5 hdf5_hl z mpi mpi_cxx curl)




##############################################################################
# UnitTest++ - http://unittest-cpp.sourceforge.net/
##############################################################################

if (NOT DEFINED UnitTest_DIR)
  message(WARNING "No UnitTest++ directory defined. To compile the unit tests add:
 -D UnitTest_DIR::FILEPATH=<UnitTest++_install_prefix>")
else()
  find_library(UNITTEST_LIBRARY 
               NAMES UnitTest++ unittest++
               HINTS ${UnitTest_LIB_DIR}
                     ${UnitTest_DIR}/lib
		     )
  if (NOT UNITTEST_LIBRARY)
    message(FATAL_ERROR "Error: could not find unittest++ library")
  endif()
  
  find_path(UNITTEST_INCLUDES "UnitTest++.h" 
            HINTS ${UnitTest_INCLUDE_DIR}
                  ${UnitTest_DIR}/include
                  ${UnitTest_DIR}/src 
                  ${UnitTest_DIR}/include/unittest++
		  )
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

##############################################################################
# Boost
##############################################################################
find_package( Boost COMPONENTS system filesystem REQUIRED)


##############################################################################
# CGNS - http://www.cgns.sourceforge.net/
##############################################################################

OPTION(ENABLE_CGNS "Build with cgnslib_2.5" OFF)

if(ENABLE_CGNS)
  find_library(CGNS_LIBRARY cgns
    HINTS ${CGNS_LIB_DIR}
    ${CGNS_DIR}/lib
    /usr/lib
    /usr/local/lib
    )
  if(NOT CGNS_LIBRARY)
    message(FATAL_ERROR "Error: could not find CGNS library\n"
      "Try defining CGNS_DIR: (-D CGNS_DIR:FILEPATH=<CGNS_install_prefix>)"
      "\nor defining CGNS_LIB_DIR")
  endif()
  
  
  find_path(CGNS_INCLUDES cgnslib.h 
    HINTS ${CGNS_INCLUDE_DIR}
    ${CGNS_DIR}/include
    /usr/include
    /usr/local/include
    )
  if(NOT CGNS_INCLUDES)
    message(FATAL_ERROR "Error: could not find CGNS includes\n"
      "Try defining CGNS_DIR: (-D CGNS_DIR:FILEPATH=<CGNS_install_prefix>)"
      "\nor defining CGNS_INCLUDE_DIR")
  endif()
  
  add_library(cgns STATIC IMPORTED)
  set_property(TARGET cgns PROPERTY IMPORTED_LOCATION ${CGNS_LIBRARY})
  set(CMAKE_CXX_FLAGS "-DENABLE_CGNS")
endif()


include_directories(${PROJECT_BINARY_DIR})

OPTION(ENABLE_DBC "Enable DBC checking" OFF)
if (ENABLE_DBC)
  message("Enabling DBC checks")
  add_definitions("-DENABLE_DBC")
endif()
  

# Enable/Disable STK_mesh Amanzi package
OPTION(ENABLE_STK_Mesh "Build Amanzi module STK_Mesh." OFF)
if (NOT ENABLE_STK_Mesh)
  message("STK_mesh compilation is disabled")
  set(MESH_LIBRARY moab_mesh)
else()
  add_subdirectory(stk_mesh)
  set(MESH_LIBRARY stk_mesh)
endif()

add_subdirectory(error_handling)
add_subdirectory(exodus)
add_subdirectory(mesh_data)
add_subdirectory(simple_mesh)
add_subdirectory(moab_mesh)
add_subdirectory(mpc)
add_subdirectory(chemistry)
add_subdirectory(transport)
add_subdirectory(boundary_conditions)
add_subdirectory(flow)
add_subdirectory(gmv)
add_subdirectory(utils)
add_subdirectory(bin)
add_subdirectory(driver)

#add_subdirectory(private)
