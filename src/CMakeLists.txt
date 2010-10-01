SET(CMAKE_CXX_COMPILER mpicxx)
SET(CMAKE_C_COMPILER mpicc)

project (STK_Mesh)
enable_testing()


cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE debug)

set(STK_Mesh_VERSION_MAJOR 0)
set(STK_Mesh_VERSION_MINOR 1)

# Curl
if(NOT DEFINED Curl_DIR)
  message(FATAL_ERROR "Error: cmake must define Curl_DIR:
(-D Curl_DIR:FILEPATH=<Curl_install_prefix>)")
  endif()

find_library(CURL_LIBRARY curl ${Curl_DIR}/lib)
if (NOT CURL_LIBRARY)
  message(FATAL_ERROR "Error: could not find Curl library")
endif()

find_path(CURL_INCLUDES curl.h ${Curl_DIR}/include/curl ${Curl_DIR}/include)
if(NOT CURL_INCLUDES)
  message(FATAL_ERROR "Error: could not find Curl includes")
endif()

add_library(curl SHARED IMPORTED)
set_property(TARGET curl PROPERTY IMPORTED_LOCATION ${CURL_LIBRARY})


# NetCDF
if(NOT DEFINED NetCDF_DIR)
  message(FATAL_ERROR "Error: cmake must define NetCDF_DIR:
 (-D NetCDF_DIR:FILEPATH=<NetCDF_install_prefix>)")
endif()

find_library(NETCDF_LIBRARY netcdf ${NetCDF_DIR}/lib)
if(NOT NETCDF_LIBRARY)
  message(FATAL_ERROR "Error: could not find NetCDF library")
endif()

find_path(NETCDF_INCLUDES netcdf.h ${NetCDF_DIR}/include)
if(NOT NETCDF_INCLUDES)
  message(FATAL_ERROR "Error: could not find NetCDF includes")
endif()

add_library(netcdf STATIC IMPORTED)
set_property(TARGET netcdf PROPERTY IMPORTED_LOCATION ${NETCDF_LIBRARY})
# target_link_libraries(netcdf curl)    See note under ExodusII.


# ExodusII
if(NOT DEFINED ExodusII_DIR)
  message(FATAL_ERROR "Error: cmake must define ExodusII_DIR:
 (-D ExodusII_DIR:FILEPATH=<ExodusII_install_prefix>)")
endif()

find_library(EXODUSII_LIBRARY exoIIv2c ${ExodusII_DIR}/lib)
if(NOT EXODUSII_LIBRARY)
  message(FATAL_ERROR "Error: could not find ExodusII library")
endif()

find_path(EXODUSII_INCLUDES exodusII.h ${ExodusII_DIR}/include)
if(NOT EXODUSII_INCLUDES)
  message(FATAL_ERROR "Error: could not find ExodusII includes")
endif()

add_library(exodusii STATIC IMPORTED)
set_property(TARGET exodusii PROPERTY IMPORTED_LOCATION ${EXODUSII_LIBRARY})

# Can't add link dependencies to imported libraries. See http://www.vtk.org/Bug/view.php?id=10395
# target_link_libraries(exodusii netcdf)

# Instead, we create a list of all the libraries required for ExodusII to link.
set(EXODUS_LIBS_ALL exodusii netcdf curl)

# UnitTest++
if (NOT DEFINED UnitTest_DIR)
  message(WARNING "No UnitTest++ directory defined. To compile the unit tests add:
 -D UnitTest_DIR::FILEPATH=<UnitTest++_install_prefix>")
else()
  find_library(UNITTEST_LIBRARY unittest++ "${UnitTest_DIR}/lib")
  if (NOT UNITTEST_LIBRARY)
    message(FATAL_ERROR "Error: could not find unittest++ library")
  endif()
  
  find_path(UNITTEST_INCLUDES "UnitTest++.h" "${UnitTest_DIR}/include" "${UnitTest_DIR}/include/unittest++")
  if (NOT UNITTEST_INCLUDES)
    message(FATAL_ERROR "Error: could not find UnitTest++ includes")
  else()
    include_directories(${UNITTEST_INCLUDES})
  endif()

  set(BUILD_TESTS 1)
  add_library(unittest++ STATIC IMPORTED)
  set_property(TARGET unittest++ PROPERTY IMPORTED_LOCATION ${UNITTEST_LIBRARY})

endif()


#Trilinos
if(NOT DEFINED Trilinos_DIR)
  message(FATAL_ERROR "Error: cmake must define Trilinos_DIR:
 (-D Trilinos_DIR:FILEPATH=<Trilinos_install_prefix>)")
endif()

find_package(Trilinos REQUIRED PATHS ${Trilinos_DIR}/include)
list(APPEND Trilinos_LIBRARIES m lapack blas)
list(APPEND Trilinos_DIRS /usr/lib)
list(APPEND Trilinos_INCLUDE_DIRS /usr/include)

include_directories(${PROJECT_BINARY_DIR})

add_subdirectory(dbc)
add_subdirectory(exodus)
add_subdirectory(mesh_data)
add_subdirectory(stk_mesh)
add_subdirectory(mpc)
add_subdirectory(boundary_conditions)
add_subdirectory(flow)
