SET(CMAKE_CXX_COMPILER mpicxx)
SET(CMAKE_C_COMPILER mpicc)

project (Amanzi)
enable_testing()


cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE debug)

set(STK_Mesh_VERSION_MAJOR 0)
set(STK_Mesh_VERSION_MINOR 1)

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


# NetCDF
find_library(NETCDF_LIBRARY netcdf
             /usr/lib
             /usr/local/lib
             ${NetCDF_DIR}/lib
             ${NetCDF_LIB_DIR})
if(NOT NETCDF_LIBRARY)
  message(FATAL_ERROR "Error: could not find NetCDF library\n"
          "Try defining NetCDF_DIR: (-D NetCDF_DIR:FILEPATH=<NetCDF_install_prefix>)"
          "\nor defining NetCDF_LIB_DIR")
endif()

find_path(NETCDF_INCLUDES netcdf.h 
          /usr/include
          /usr/local/include
          ${NetCDF_DIR}/include
          ${NetCDF_INCLUDE_DIR})
if(NOT NETCDF_INCLUDES)
  message(FATAL_ERROR "Error: could not find NetCDF includes\n"
          "Try defining NetCDF_DIR: (-D NetCDF_DIR:FILEPATH=<NetCDF_install_prefix>)"
          "\nor defining NetCDF_INCLUDE_DIR")
endif()

add_library(netcdf STATIC IMPORTED)
set_property(TARGET netcdf PROPERTY IMPORTED_LOCATION ${NETCDF_LIBRARY})
# target_link_libraries(netcdf curl)    See note under ExodusII.


# ExodusII
find_library(EXODUSII_LIBRARY exoIIv2c ${ExodusII_DIR}/lib)
if(NOT EXODUSII_LIBRARY)
  message(FATAL_ERROR "Error: could not find ExodusII library\n"
          "Try defining ExodusII_DIR: (-D ExodusII_DIR:FILEPATH=<ExodusII_install_prefix>)\n"
          "Or ExodusII_LIB_DIR")
endif()

find_path(EXODUSII_INCLUDES exodusII.h 
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
add_subdirectory(mpc)
add_subdirectory(chemistry)
add_subdirectory(transport)
add_subdirectory(boundary_conditions)
add_subdirectory(flow)
add_subdirectory(bin)
