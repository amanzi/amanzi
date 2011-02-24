# -*- mode: cmake -*-
    
project (AMANZI)
enable_testing()

set(ARCHOS ${CMAKE_SYSTEM_PROCESSOR}_${CMAKE_SYSTEM_NAME})

cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE debug)

set(AMANZI_Mesh_VERSION_MAJOR 0)
set(AMANZI_Mesh_VERSION_MINOR 1)

if(PREFER_STATIC_LIBRARIES)
  # Prefer static libraries, but don't require that everything must be static. 
  # This appears to be necessary on Franklin at NERSC right now.  --RTM
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a .lib)
endif(PREFER_STATIC_LIBRARIES)

if(BUILD_STATIC_EXECUTABLES)
    set(CMAKE_EXE_LINKER_FLAGS -static)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)       # remove -Wl,-Bdynamic
    set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
    set(CMAKE_SHARED_LIBRARY_C_FLAGS)         # remove -fPIC
    set(CMAKE_SHARED_LIBRARY_CXX_FLAGS)
    set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)    # remove -rdynamic
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
endif(BUILD_STATIC_EXECUTABLES)



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

#target_link_libraries(netcdf)   # See note under ExodusII.
set(NETCDF_LIBS_ALL netcdf hdf5 hdf5_hl z)

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

set(HDF5_LIBS_ALL hdf5 hdf5_hl z)



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
set(EXODUS_LIBS_ALL exodusii netcdf hdf5 hdf5_hl z ) 




##############################################################################
# MOAB - svn co https://svn.mcs.anl.gov/repos/ITAPS/MOAB/trunk MOAB
##############################################################################

OPTION(ENABLE_MOAB_Mesh "Build with MOAB" OFF)
if (ENABLE_MOAB_Mesh)
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
  set(MOAB_LIBS_ALL moab netcdf_c++ netcdf hdf5 hdf5_hl z )

  endif(ENABLE_MOAB_Mesh)



##############################################################################
# MSTK - https://software.lanl.gov/MeshTools/trac/raw-attachment/wiki/WikiStart/mstk-1.80.tar.gz
##############################################################################

OPTION(ENABLE_MSTK_Mesh "Build Amanzi module MSTK_Mesh." OFF)
if (ENABLE_MSTK_Mesh)
  find_library(MSTK_LIBRARY mstk
             HINTS ${MSTK_LIB_DIR}
                   ${MSTK_DIR}/lib
                   ${MSTK_DIR}/lib/${ARCHOS}
                   /usr/lib
                   /usr/local/lib
		   )
  if(NOT MSTK_LIBRARY)
    message(FATAL_ERROR "Error: could not find MSTK library\n"
          "Try defining MSTK_DIR: (-D MSTK_DIR:FILEPATH=<MSTK_install_prefix>)"
          "\nor defining MSTK_LIB_DIR")
  endif()

  find_path(MSTK_INCLUDES MSTK.h 
          HINTS ${MSTK_INCLUDE_DIR}
                ${MSTK_DIR}/include
                /usr/include
                /usr/local/include
		)
  if(NOT MSTK_INCLUDES)
    message(FATAL_ERROR "Error: could not find MSTK includes\n"
          "Try defining MSTK_DIR: (-D MSTK_DIR:FILEPATH=<MSTK_install_prefix>)"
          "\nor defining MSTK_INCLUDE_DIR")
  endif()

  add_library(mstk STATIC IMPORTED)
  set_property(TARGET mstk PROPERTY IMPORTED_LOCATION ${MSTK_LIBRARY})

# create a list of all the libraries required for MSTK to link.
  set(MSTK_LIBS_ALL mstk exodusii netcdf metis)


##############################################################################
# METIS - http://glaros.dtc.umn.edu/gkhome/metis/metis/download
##############################################################################

  find_library(METIS_LIBRARY metis
	     HINTS ${METIS_LIB_DIR}
                   ${METIS_DIR}
	           ${METIS_DIR}/lib
		   )
  if(NOT METIS_LIBRARY)
    message(FATAL_ERROR "Error: could not find METIS library\n"
          "Try defining METIS_DIR: (-D METIS_DIR:FILEPATH=<Metis_install_prefix>)\n"
          "Or METIS_LIB_DIR")
  endif()

  find_path(METIS_INCLUDES metis.h 
          HINTS ${METIS_INCLUDE_DIR}
                ${METIS_DIR}
	        ${METIS_DIR}/include
		${METIS_DIR}/Lib
		)
  if(NOT METIS_INCLUDES)
    message(FATAL_ERROR "Error: could not find METIS includes\n"
          "Try defining METIS_DIR: (-D METIS_DIR:FILEPATH=<METIS_install_prefix>)\n"
          "Or METIS_INCLUDE_DIR")
  endif()

  add_library(metis STATIC IMPORTED)
  set_property(TARGET metis PROPERTY IMPORTED_LOCATION ${METIS_LIBRARY})

endif(ENABLE_MSTK_Mesh)


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

find_package(Trilinos REQUIRED PATHS ${Trilinos_DIR}/include)
if(NOT Trilinos_FOUND)
  message(FATAL_ERROR "Error: could not find Trilinos\n"
          "Try defining Trilinos_DIR: (-D Trilinos_DIR:FILEPATH=<Trilinos_install_prefix>)")
endif()
list(APPEND Trilinos_DIRS /usr/lib)
list(APPEND Trilinos_INCLUDE_DIRS /usr/include)

##############################################################################
# Boost
##############################################################################
find_package( Boost COMPONENTS system filesystem program_options REQUIRED)


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
  add_definitions("-DENABLE_CGNS")
else()
  message(WARNING "CGNS not enabled. Some modules and tests will not be compiled. Use option"
    "-D ENABLE_CGNS::BOOL=True to enable.")
endif()


include_directories(${PROJECT_BINARY_DIR})

OPTION(ENABLE_DBC "Enable DBC checking" ON)
if (ENABLE_DBC)
  message("Enabling DBC checks")
  add_definitions("-DENABLE_DBC")
endif()
  
##############################################################################
# STK_mesh Amanzi package
##############################################################################
OPTION(ENABLE_STK_Mesh "Build Amanzi module STK_Mesh." OFF)
if (NOT ENABLE_STK_Mesh)
  message("STK_mesh compilation is disabled")
else()
  add_subdirectory(stk_mesh)
  set(MESH_LIBRARY stk_mesh)
  add_definitions("-DENABLE_STK")
endif()

if (ENABLE_MSTK_Mesh) 
  add_subdirectory(mstk_mesh)
  set(MESH_LIBRARY mstk_mesh)
endif()

if (ENABLE_MOAB_Mesh)
  add_subdirectory(moab_mesh)
endif()

add_subdirectory(error_handling)
add_subdirectory(output)
add_subdirectory(exodus)
add_subdirectory(mesh)
add_subdirectory(mesh_data)
add_subdirectory(simple_mesh)
add_subdirectory(mpc)
add_subdirectory(chemistry)
add_subdirectory(transport)
add_subdirectory(boundary_conditions)
add_subdirectory(flow)
add_subdirectory(bin)
add_subdirectory(driver)


message(STATUS "******************************************************************************")
message(STATUS "")
message(STATUS "${PROJECT_NAME} Configuration")
message(STATUS "Version ${AMANZI_VERSION_MAJOR}.${AMANZI_VERSION_MINOR}")
message(STATUS "")
message(STATUS "System Information")
message(STATUS "System ${CMAKE_SYSTEM}")
message(STATUS "System Name ${CMAKE_SYSTEM_NAME}")
message(STATUS "System Version ${CMAKE_SYSTEM_VERSION}")
if ( APPLE )
    message(STATUS "Mac OSX system")
endif(APPLE)    
if ( WIN32 )
    message(STATUS "Windows system")
endif(WIN32)    
if ( UNIX )
    message(STATUS "UNIX system")
endif(UNIX)   
message(STATUS "")
message(STATUS "Build type ${CMAKE_BUILD_TYPE}")
message(STATUS "")
message(STATUS "Compilers")
message(STATUS "C COMPILER ${CMAKE_C_COMPILER}")
message(STATUS "CXX COMPILER ${CMAKE_CXX_COMPILER}")
message(STATUS "")
message(STATUS "******************************************************************************")
