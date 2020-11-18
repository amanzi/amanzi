#  -*- mode: cmake -*-

#
# Amanzi SuperBuild
#       Amanzi TPL Build Targets

#
# Set minimum required version of CMake
# NB: Must go before project() call.
# KL: Version used to build Amanzi must be the same!
cmake_minimum_required(VERSION 3.11.4)

#
# Language standard
# NB: Seems like a good idea to have this before project() call.
#
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ############################################################################ #
# Define project name (SuperBuild_SOURCE_DIR and SuperBuild_BINARY_DIR)
# ############################################################################ #
message(STATUS "\n------------------------------------------------")
message(STATUS "\n-- CMake: Configuring Amanzi TPLs build/install.\n--")
message(STATUS "---------------------------------------------")
project(SuperBuild C CXX Fortran)

# ############################################################################ #
# Shared Libraries Definitions
# ############################################################################ #
if (BUILD_SHARED_LIBS)
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_SKIP_INSTALL_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH "${TPL_INSTALL_PREFIX}/lib")
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE) 
  set(AMANZI_CMAKE_CACHE_ARGS  "-DCMAKE_INSTALL_RPATH:PATH=${CMAKE_INSTALL_RPATH}")
  list(APPEND AMANZI_CMAKE_CACHE_ARGS  "-DCMAKE_INSTALL_NAME_DIR:PATH=${CMAKE_INSTALL_RPATH}")
else ()
  SET(CMAKE_SKIP_INSTALL_RPATH ON CACHE BOOL "Turn off for static install." FORCE)
  SET(CMAKE_SKIP_RPATH ON CACHE BOOL "Turn off for static install." FORCE)
endif()

if (BUILD_SHARED_LIBS)
  set(_shared_flags "-fPIC")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${_shared_flags}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_shared_flags}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${_shared_flags}")

  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${_shared_flags}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${_shared_flags}")
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${_shared_flags}")

  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${_shared_flags}")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${_shared_flags}")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${_shared_flags}")
endif()

# ############################################################################ #
# Static Libraries Definitions
# ############################################################################ #
set(BUILD_STATIC_LIBS TRUE)
if (BUILD_SHARED_LIBS)
  set(BUILD_STATIC_LIBS FALSE)
endif()


# ############################################################################ #
# Compiler Definition Warnings
# ############################################################################ #
if (NOT CMAKE_C_COMPILER)
  message(WARNING "C compiler not specified. CMake will guess!")
endif()

if (NOT CMAKE_CXX_COMPILER)
  message(WARNING "C++ compiler not specified. CMake will guess!")
endif()

if (NOT CMAKE_Fortran_COMPILER)
  message(WARNING "Fortran compiler not specified. CMake will guess!")
endif()

#if (NOT SPACK_BINARY)
#  message(WARNING "Spack binary not defined. Disabling spack builds.")
#endif()

# ############################################################################ #
# Define Amanzi directories 
# ############################################################################ #
file(TO_CMAKE_PATH ${SuperBuild_SOURCE_DIR}/../.. Amanzi_SOURCE_DIR)

# ############################################################################ #
# Update the CMake module path
# ############################################################################ #
set(Amanzi_MODULE_PATH ${Amanzi_SOURCE_DIR}/tools/cmake
                       ${Amanzi_SOURCE_DIR}/tools/cmake/Utils
                       ${Amanzi_SOURCE_DIR}/tools/cmake/Modules)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} 
	              ${SuperBuild_SOURCE_DIR}/cmake
	              ${Amanzi_MODULE_PATH})

# ############################################################################ #
# Location of the Build_* files and template files
# ############################################################################ #
set(SuperBuild_BUILD_FILES_DIR    "${SuperBuild_SOURCE_DIR}/include")
set(SuperBuild_TEMPLATE_FILES_DIR "${SuperBuild_SOURCE_DIR}/templates")

# ############################################################################ #
# Compiler checks
# ############################################################################ #
include(FortranCInterface)
FortranCInterface_VERIFY()

include(DefineCompilerVersion)
define_compiler_version()

include(CheckMPISourceCompiles)
check_mpi_source_compiles(MPI_WRAPPERS_IN_USE)
if (NOT MPI_WRAPPERS_IN_USE)
  # MPI_WRAPPERS_IN_USE is FALSE if any of the compilers fail to compile an MPI app
  message(FATAL_ERROR "Missing a compiler type for MPI applications -- exiting.")
endif()

# Include the TPL version information
include(${SuperBuild_SOURCE_DIR}/TPLVersions.cmake)

# prepare the version numbers file
set(TPL_VERSIONS_INCLUDE_FILE ${CMAKE_BINARY_DIR}/tpl_versions.h)
file(WRITE ${TPL_VERSIONS_INCLUDE_FILE} "")
install(FILES ${CMAKE_BINARY_DIR}/tpl_versions.h DESTINATION include)
amanzi_tpl_version_write(FILENAME ${TPL_VERSIONS_INCLUDE_FILE}
  PREFIX AMANZI
  VERSION ${AMANZI_TPLS_VERSION_MAJOR} ${AMANZI_TPLS_VERSION_MINOR} ${AMANZI_TPLS_VERSION_PATCH})


# ############################################################################ #
# Build options
# ############################################################################ #
set(SuperBuild_BUILD_TYPE_DFLT "Release")
if (NOT CMAKE_BUILD_TYPE) 
  set(CMAKE_BUILD_TYPE ${SuperBuild_BUILD_TYPE_DFLT})
endif()
message(STATUS "Setting build type -- ${CMAKE_BUILD_TYPE}")

# Download directory TPL_DOWNLOAD_DIR
set(TPL_DOWNLOAD_DIR_DFLT ${SuperBuild_BINARY_DIR}/Downloads)
if (NOT TPL_DOWNLOAD_DIR)
  set(TPL_DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR_DFLT})
endif()
message(STATUS "Download TPL files in ${TPL_DOWNLOAD_DIR}") 

# Disable external web downloads
option(DISABLE_EXTERNAL_DOWNLOAD "Flag to disable external web downloads" FALSE)
if (DISABLE_EXTERNAL_DOWNLOAD)
  message(STATUS "Will not download files from the web, instead will search "
                 "for archive files in ${TPL_DOWNLOAD_DIR}")
else()
  include(CheckDownload)
  check_download(TEST_FILE MD5SUMS TEST_URL ${AMANZI_TPLS_DOWNLOAD_URL})
endif()

# TPL install location
set(TPL_INSTALL_PREFIX_DFLT "${SuperBuild_BINARY_DIR}/external-projects")
if (NOT TPL_INSTALL_PREFIX)
  set(TPL_INSTALL_PREFIX ${TPL_INSTALL_PREFIX_DFLT})
endif()
message(STATUS "Install TPLs in ${TPL_INSTALL_PREFIX}")

# also set the CMAKE_INSTALL_PREFIX accordingly
set(CMAKE_INSTALL_PREFIX ${TPL_INSTALL_PREFIX})

# Amanzi structured AMR mesh API and tools
option(ENABLE_Structured "Enable Amanzi structured mesh capability" TRUE)
option(ENABLE_CCSE_TOOLS "Enable building CCSE post processing tools" FALSE)

# Amanzi unstructured mesh API and frameworks
option(ENABLE_Unstructured "Enable Amanzi unstructured mesh capability" TRUE)
option(ENABLE_STK_Mesh "Enable the Unstructured STK Mesh Toolkit" FALSE)
option(ENABLE_MSTK_Mesh "Enable the Unstructured MSTK Mesh Toolkit" TRUE)
option(ENABLE_MOAB_Mesh "Enable the Unstructured MOAB Mesh Toolkit" FALSE)

# Amanzi output formats
option(ENABLE_Silo "Enable the Silo output specification for visualization" FALSE)

# OpenMP Search - ENABLED ONLY WHEN STRUCTURED IS ON BECAUSE
# WE SUSPECT THAT IT MAY NOT BE WORKING WELL WITH HYPRE ON
# THE UNSTRUCTURED SIDE
if (ENABLE_Structured)
  # option(ENABLE_OpenMP "Add OpenMP directives to the build" TRUE)
  set(ENABLE_OpenMP FALSE)
  if (ENABLE_OpenMP)
    find_package(OpenMP)
    find_package(OpenMP_Fortran)
    if (OPENMP_FOUND)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
      set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${OpenMP_LINK_FLAGS}")
    else()  
      message(SEND_ERROR "Failed to determine the OpenMP compiler flags. Will deactivate.")
      set(ENABLE_OpenMP FALSE CACHE BOOL "Enable OpenMP compile directives" FORCE)
    endif()
  endif()
endif()

# If unstructured is not enabled, there is no reason to have these mesh frameworks
if (NOT ENABLE_Unstructured)
  set(ENABLE_STK_Mesh FALSE)
  set(ENABLE_MSTK_Mesh FALSE)
  set(ENABLE_MOAB_Mesh FALSE)
endif()

# Amanzi Test Suite
option(ENABLE_TESTS "Build the Amanzi Test Suite" TRUE)

# xSDK Libraries
option(ENABLE_XSDK "Use xSDK Libraries" FALSE)

# Use Spack package manager
option(ENABLE_SPACK "Build with the Spack package manager where appropriate" FALSE)
option(BUILD_SPACK "Clone spack from git and treat as a TPL that needs to be built" FALSE)


# ############################################################################ #
# MPI Configuration 
# ############################################################################ #

# --- The FindMPI package for CMake for versions 2.8.6> use compiler wrappers
#     to determine include paths, link flags and libraries for each language.
if (MPI_WRAPPERS_IN_USE)
  set(MPI_C_COMPILER ${CMAKE_C_COMPILER})
  set(MPI_CXX_COMPILER ${CMAKE_CXX_COMPILER})
  set(MPI_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})

else()
  message(STATUS "Searching for MPI compiler wrappers to define"
                 " compile and link flags")

  find_program(MPI_C_COMPILER
               NAME mpicc
               HINTS ${MPI_INSTALL_PREFIX} ENV MPIROOT ENV MPI_ROOT ENV MPI_HOME ENV MPIHOME ENV MPICH_ROOT 
               PATH_SUFFIXES bin Bin
               )

  find_program(MPI_CXX_COMPILER
               NAMES mpicxx mpic++ mpiCC
               HINTS ${MPI_INSTALL_PREFIX} ENV MPIROOT ENV MPI_ROOT ENV MPI_HOME ENV MPIHOME ENV MPICH_ROOT 
               PATH_SUFFIXES bin Bin
               )

  find_program(MPI_Fortran_COMPILER
               NAMES mpif90
               HINTS ${MPI_INSTALL_PREFIX} ENV MPIROOT ENV MPI_ROOT ENV MPI_HOME ENV MPIHOME ENV MPICH_ROOT 
               PATH_SUFFIXES bin Bin
               )
endif()

# --- This package find call will set
#     MPI_<lang>_FOUND
#     MPI_<lang>_COMPILE_FLAGS
#     MPI_<lang>_INCLUDE_PATH
#     MPI_<lang>_LINK_FLAGS
#     MPI_<lang>_LIBRARIES
#
#     MPIEXEC
#     MPIEXEC_NUMPROC_FLAG
#     MPIEXEC_PREFLAGS
#     MPIEXEC_POSTFLAGS
find_package(MPI)


# The FindMPI package uses FindPackageHandleStandardArgs
# to MPI_<lang>_FOUND. For consistency, I set MPI_FORTRAN_FOUND
# to MPI_Fortran_FOUND
set(MPI_Fortran_FOUND ${MPI_FORTRAN_FOUND})

# Warn the user if MPI information is not found
if (NOT MPI_C_FOUND)
  message(WARNING "Failed to locate MPI C include and library files")
endif()

if (NOT MPI_CXX_FOUND)
  message(WARNING "Failed to locate MPI C++ include and library files")
endif()

if (NOT MPI_Fortran_FOUND)
  message(WARNING "Failed to locate MPI Fortran include and library files")
endif()

# Fail configuration if any MPI_*_FOUND == FALSE
if (MPI_C_FOUND AND MPI_CXX_FOUND AND MPI_Fortran_FOUND)
  message(STATUS "Found MPI: MPIEXEC=${MPIEXEC}")
  set(mpi_var_list COMPILER COMPILE_FLAGS INCLUDE_PATH LINK_FLAGS LIBRARIES)
  foreach(lang C CXX Fortran)
    foreach(var ${mpi_var_list})
      set(var MPI_${lang}_${var})
      if ("${var}")
	message(STATUS "      ${var}: ${${var}}")
      endif()
    endforeach()
  endforeach()  
else()
  message(FATAL_ERROR "MPI was not found.")
endif()

# Add the include paths, link flags and libraries if using a 
# user requested MPI but not the wrappers from this MPI installation.
if (NOT MPI_WRAPPERS_IN_USE)
  foreach (lang C CXX Fortran)
    # Build MPI compile flags
    set(add_mpi_flags_list)
    if (MPI_${lang}_COMPILE_FLAGS)
      foreach (flag ${MPI_${lang}_COMPILE_FLAGS} )
	list(APPEND add_mpi_flags_list ${flag})
      endforeach()
    endif()

    if (MPI_${lang}_INCLUDE_PATH)
      foreach (dir ${MPI_${lang}_INCLUDE_PATH} )
	list(APPEND add_mpi_flags_list "-I${dir}")
      endforeach()
    endif() 

    if (MPI_${lang}_LINK_FLAGS)
      foreach (flag ${MPI_${lang}_LINK_FLAGS} )
	list(APPEND add_mpi_flags_list ${flag})
      endforeach()
    endif()

    list(REMOVE_DUPLICATES add_mpi_flags_list)
    list(LENGTH add_mpi_flags_list num_flags)

    if ("${num_flags}" GREATER "0")
      include(BuildWhitespaceString)
      if (CMAKE_${lang}_FLAGS)
	build_whitespace_string(new_compile_flags ${CMAKE_${lang}_FLAGS} ${add_mpi_flags_list})
      else()	
	build_whitespace_string(new_compile_flags ${add_mpi_flags_list})
      endif()
      set(CMAKE_${lang}_FLAGS "${new_compile_flags}")
      message(STATUS "Update CMAKE_${lang}_FLAGS: ${CMAKE_${lang}_FLAGS}")
    endif()  
  endforeach()
endif()

# --- Amanzi uses MPI_EXEC* not MPIEXEC* variables. This allows the user to 
#     override the find package results.

# - MPI execute binary
if (NOT MPI_EXEC)
  if (MPIEXEC)
    set(MPI_EXEC "${MPIEXEC}" CACHE STRING "Select MPI executable from FindMPI results")
  else()
    find_program(MPI_EXEC
                 NAMES mpirun mpiexec aprun openmpirun
                 HINTS ENV MPI_ROOT ENV MPIROOT ENV MPI_PREFIX ENV MPI_HOME ENV MPIHOME
                )
  endif()
endif()

# - Number of MPI ranks flag
set(MPI_EXEC_NUMPROCS_FLAG_DFLT -n)
if (NOT MPI_EXEC_NUMPROCS_FLAG)
  if (MPIEXEC_NUMPROC_FLAG)
    set(MPI_EXEC_NUMPROCS_FLAG "${MPIEXEC_NUMPROC_FLAG}" CACHE STRING "Set MPI number of procs flag from FindMPI")
  else()
    set(MPI_EXEC_NUMPROCS_FLAG ${MPI_EXEC_NUMPROCS_FLAG_DFLT})
  endif()
endif()

# - Maximum number of processors. This is a limit for the test suite
#   Some tests require too many processors and it increases the execution time
#   considerably. 
set(MPI_EXEC_MAX_NUMPROCS_DFLT 8)
if (NOT MPI_EXEC_MAX_NUMPROCS)
  include(ProcessorCount)
  ProcessorCount(proc_count)
  if (NOT proc_count EQUAL 0)
    math(EXPR MPI_EXEC_MAX_NUMPROCS "${proc_count} * 2") 
    message(STATUS "Detected ${proc_count} processors and will set maximum to ${MPI_EXEC_MAX_NUMPROCS}")
  else()
    set(MPI_EXEC_MAX_NUMPROCS ${MPI_EXEC_MAX_NUMPROCS_DFLT})
  endif()
endif()  

# - Set the pre and post flags
#   Usage:
#   ${MPI_EXEC} ${MPI_EXEC_NUMPROCS_FLAG} PROCS ${MPI_EXEC_PREFLAGS} EXECUTABLE ${MPI_EXEC_POSTFLAGS}
if (NOT MPI_EXEC_PREFLAGS)
  if (MPIEXEC_PREFLAGS)
    set(MPI_EXEC_PREFLAGS "${MPIEXEC_PRFLAGS}" CACHE STRING "Set MPI execute pre flags")
  endif()
endif()

if (NOT MPI_EXEC_POSTFLAGS)
  if (MPIEXEC_POSTFLAGS)
    set(MPI_EXEC_POSTFLAGS "${MPIEXEC_PRFLAGS}" CACHE STRING "Set MPI execute post flags")
  endif()
endif()


# ############################################################################ #
# LAPACK/BLAS Configuration 
# ############################################################################ #

# Search for BLAS (Basic Linear Algebra Subprograms) http://www.netlib.org/blas
find_package(BLAS REQUIRED)
message(STATUS ">>>  BLAS_LIBRARIES = ${BLAS_LIBRARIES}") 

# Search for LAPACK (Linear Algebra PACKage) http://www/netlib.org/lapack
find_package(LAPACK REQUIRED)
message(STATUS ">>>  LAPACK_LIBRARIES = ${LAPACK_LIBRARIES}") 


# ############################################################################ #
# Set common build compiler flags, build types and directories
# ############################################################################ #

# Common compiler flags to ensure a uniform build in projects that do not use CMake
include(BuildWhitespaceString)
string(TOUPPER "${CMAKE_BUILD_TYPE}" build_type_uc)
set(build_c_flags ${CMAKE_C_FLAGS_${build_type_uc}})
set(build_cxx_flags ${CMAKE_CXX_FLAGS_${build_type_uc}})
set(build_fortran_flags ${CMAKE_Fortran_FLAGS_${build_type_uc}})
message(STATUS "Common compile flags for build type ${CMAKE_BUILD_TYPE}")
build_whitespace_string(Amanzi_COMMON_CFLAGS ${CMAKE_C_FLAGS} ${build_c_flags}) 
build_whitespace_string(Amanzi_COMMON_CXXFLAGS ${CMAKE_CXX_FLAGS} ${build_cxx_flags}) 
build_whitespace_string(Amanzi_COMMON_FCFLAGS ${CMAKE_Fortran_FLAGS} ${build_fortran_flags}) 

# Common link flags for projects that do not use CMake
set(Amanzi_COMMON_LDFLAGS)
if (CMAKE_EXE_LINK_FLAGS)
  build_whitespace_string(Amanzi_COMMON_LDFLAGS ${CMAKE_EXE_LINK_FLAGS})
endif()

# Common CMake controls
list(APPEND AMANZI_CMAKE_CACHE_ARGS "-DCMAKE_PREFIX_PATH:PATH=${TPL_INSTALL_PREFIX}")

message(STATUS "\tcommon C flags\t\t${Amanzi_COMMON_CFLAGS}")
message(STATUS "\tcommon C++ flags\t${Amanzi_COMMON_CXXFLAGS}")
message(STATUS "\tcommon Fortran flags\t${Amanzi_COMMON_FCFLAGS}")
message(STATUS "\tcommon Linker flags\t${Amanzi_COMMON_LDFLAGS}")
message(STATUS "\tcommon arguments\t${AMANZI_CMAKE_CACHE_ARGS}")


# ############################################################################ #
# Begin TPL builds 
# ############################################################################ #

# --- Add .NOTPARALLEL target in order to force the top level make to 
#     execute serially

add_custom_target(.NOTPARALLEL)

# --- Initalize and set the TPL build files

# Create configure, build, install and test targets for each TPL
set_property(DIRECTORY PROPERTY
             EP_STEP_TARGETS download patch configure build install test)

# Include the CMake module ExternalProject
include(ExternalProject)

# Include the macro that defines uniform build, install and Add_ExternalProject args
include(DefineExternalProjectArgs)

# Set of CMake configure arguments to pass to Amanzi once the TPLs are built
set(Amanzi_TPL_CMAKE_ARGS)

# --- Begin the build definitions
#if (BUILD_SPACK)
#  include(${SuperBuild_BUILD_FILES_DIR}/Build_Spack.cmake)
#endif()


if (ENABLE_XSDK)    
  include(${SuperBuild_BUILD_FILES_DIR}/Build_xSDK.cmake)
  include(BuildLibraryName)
  include(FindSpackInstallLocation)
endif()

# ZLIB
if (NOT ENABLE_XSDK)
  include(${SuperBuild_BUILD_FILES_DIR}/Build_zlib.cmake)
else()
  find_spack_install_location(zlib ZLIB_INSTALL_PREFIX)
  message(STATUS "ZLIB_INSTALL_PREFIX = ${ZLIB_INSTALL_PREFIX}")
  build_library_name(z ZLIB_LIBRARIES APPEND_PATH ${ZLIB_INSTALL_PREFIX}/lib)
  set(ZLIB_INCLUDE_DIRS ${ZLIB_INSTALL_PREFIX}/include)
  set(ZLIB_DIR ${ZLIB_INSTALL_PREFIX})
endif()

# UnitTest
append_set(Amanzi_TPL_CMAKE_ARGS -DENABLE_TEST:BOOL=${ENABLE_TESTS})
if (ENABLE_TESTS)
  include(${SuperBuild_BUILD_FILES_DIR}/Build_UnitTest.cmake)
  append_set(Amanzi_TPL_CMAKE_ARGS 
             -DUnitTest_DIR:FILEPATH=${TPL_INSTALL_PREFIX})
endif()

# CCSE
set(CCSE_BL_SPACEDIM_DFLT 2)
if (NOT CCSE_BL_SPACEDIM )
  set(CCSE_BL_SPACEDIM ${CCSE_BL_SPACEDIM_DFLT})
endif()
if (ENABLE_Structured)
  include(${SuperBuild_BUILD_FILES_DIR}/Build_CCSE.cmake) 
  append_set(Amanzi_TPL_CMAKE_ARGS
             -DCCSE_DIR:FILEPATH=${TPL_INSTALL_PREFIX}
             -DENABLE_MPI:BOOL=TRUE
             -DENABLE_OpenMP:BOOL=${ENABLE_OpenMP}
             -DAMANZI_PRECISION:STRING=DOUBLE
             -DAMANZI_SPACEDIM:INT=${CCSE_BL_SPACEDIM})
endif()

# Boost
if (NOT ENABLE_XSDK)
    include(${SuperBuild_BUILD_FILES_DIR}/Build_Boost.cmake) 
    append_set(Amanzi_TPL_CMAKE_ARGS
        -DBOOST_ROOT:FILEPATH=${TPL_INSTALL_PREFIX}
        -DBOOST_INCLUDEDIR=${TPL_INSTALL_PREFIX}/include/boost
        -DBoost_NO_SYSTEM_PATHS:BOOL=TRUE)
else()
    find_spack_install_location(boost BOOST_INSTALL_PREFIX)
    message(STATUS "BOOST_INSTALL_PREFIX = ${BOOST_INSTALL_PREFIX}")
    set(BOOST_ROOT ${BOOST_INSTALL_PREFIX})
    set(BOOST_INCLUDEDIR ${BOOST_INSTALL_PREFIX}/include/boost)
    append_set(Amanzi_TPL_CMAKE_ARGS
        -DBOOST_ROOT:FILEPATH=${BOOST_INSTALL_PREFIX}
        -DBOOST_INCLUDEDIR=${BOOST_INSTALL_PREFIX}/include/boost
        -DBoost_USE_STATIC_LIBS:BOOL=TRUE
        -DBoost_NO_SYSTEM_PATHS:BOOL=TRUE)
endif()

# HDF5  
if (NOT ENABLE_XSDK)
  include(${SuperBuild_BUILD_FILES_DIR}/Build_HDF5.cmake)
else()
  amanzi_tpl_version_write(FILENAME ${TPL_VERSIONS_INCLUDE_FILE}
                           PREFIX HDF5
                           VERSION 1 10 0)
  find_spack_install_location(hdf5 HDF5_INSTALL_PREFIX)
  message(STATUS "HDF5_INSTALL_PREFIX = ${HDF5_INSTALL_PREFIX}")
  build_library_name(hdf5 HDF5_C_LIBRARY APPEND_PATH ${HDF5_INSTALL_PREFIX}/lib)
  build_library_name(hdf5_hl HDF5_HL_LIBRARY APPEND_PATH ${HDF5_INSTALL_PREFIX}/lib)
  set(HDF5_LIBRARIES ${HDF5_HL_LIBRARY} ${HDF5_C_LIBRARY} ${ZLIB_LIBRARIES} m)
  set(HDF5_INCLUDE_DIRS ${HDF5_INSTALL_PREFIX}/include ${ZLIB_INCLUDE_DIRS})
  set(HDF5_DIR ${HDF5_INSTALL_PREFIX})
endif()
append_set(Amanzi_TPL_CMAKE_ARGS
           -DHDF5_ROOT:FILEPATH=${HDF5_DIR})

# ASCEM-IO  
include(${SuperBuild_BUILD_FILES_DIR}/Build_ASCEMIO.cmake) 
append_set(Amanzi_TPL_CMAKE_ARGS
           -DASCEMIO_DIR:FILEPATH=${TPL_INSTALL_PREFIX})

# XERCES
include(${SuperBuild_BUILD_FILES_DIR}/Build_XERCES.cmake) 
append_set(Amanzi_TPL_CMAKE_ARGS
          -DXERCES_DIR:FILEPATH=${XERCES_DIR})

# NetCDF
if (NOT ENABLE_XSDK)
    include(${SuperBuild_BUILD_FILES_DIR}/Build_NetCDF.cmake) 
else()
    find_spack_install_location(netcdf NETCDF_INSTALL_PREFIX)
    message(STATUS "NETCDF_INSTALL_PREFIX = ${NETCDF_INSTALL_PREFIX}")
    build_library_name(netcdf NetCDF_C_LIBRARY APPEND_PATH ${NETCDF_INSTALL_PREFIX}/lib)
    build_library_name(netcdf_c++ NetCDF_CXX_LIBRARY APPEND_PATH ${NETCDF_INSTALL_PREFIX}/lib)
    set(NetCDF_INCLUDE_DIRS ${NETCDF_INSTALL_PREFIX}/include)
    set(NetCDF_DIR ${NETCDF_INSTALL_PREFIX})
    set(NetCDF_C_LIBRARIES ${NetCDF_C_LIBRARY})
    if (ENABLE_NetCDF4)
      list(APPEND NetCDF_C_LIBRARIES ${HDF5_LIBRARIES})
      list(APPEND NetCDF_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
      list(REMOVE_DUPLICATES NetCDF_INCLUDE_DIRS)
    endif()
endif()
append_set(Amanzi_TPL_CMAKE_ARGS
          -DNetCDF_DIR:FILEPATH=${NETCDF_DIR})

# NetCDF Fortran
# As of 10/13/2017, NetCDF Fortran is not part of the xSDK suite
include(${SuperBuild_BUILD_FILES_DIR}/Build_NetCDF_Fortran.cmake)
append_set(Amanzi_TPL_CMAKE_ARGS
          -DNetCDF_FORTRAN_DIR:FILEPATH=${NetCDF_FORTRAN_DIR})

# MOAB
if (ENABLE_MOAB_Mesh)
  include(${SuperBuild_BUILD_FILES_DIR}/Build_MOAB.cmake)
  append_set(Amanzi_TPL_CMAKE_ARGS
             -DENABLE_MOAB_Mesh:BOOL=${ENABLE_MOAB_Mesh}
             -DMOAB_DIR:PATH=${TPL_INSTALL_PREFIX})
endif()
 
# METIS -- needs to go before PETSc because PETSc needs it.
if (NOT ENABLE_XSDK)
    include(${SuperBuild_BUILD_FILES_DIR}/Build_METIS.cmake) 
else()
    find_spack_install_location(metis METIS_INSTALL_PREFIX)
    message(STATUS "METIS_INSTALL_PREFIX = ${METIS_INSTALL_PREFIX}")
    build_library_name(metis METIS_LIBRARIES SHARED APPEND_PATH ${METIS_INSTALL_PREFIX}/lib)
    set(METIS_INCLUDE_DIR ${METIS_INSTALL_PREFIX}/include)
    set(METIS_DIR ${METIS_INSTALL_PREFIX})
endif()
append_set(Amanzi_TPL_CMAKE_ARGS 
           -DMETIS_DIR:FILEPATH=${METIS_DIR})

# ParMetis -- needs to go before PETSc because PETSc needs it.
if (NOT ENABLE_XSDK)
    include(${SuperBuild_BUILD_FILES_DIR}/Build_ParMetis.cmake)
else()
    find_spack_install_location(parmetis PARMETIS_INSTALL_PREFIX)
    message(STATUS "PARMETIS_INSTALL_PREFIX = ${PARMETIS_INSTALL_PREFIX}")
    build_library_name(z ParMetis_LIBRARIES APPEND_PATH ${PARMETIS_INSTALL_PREFIX}/lib)
    set(ParMetis_INCLUDE_DIRS ${PARMETIS_INSTALL_PREFIX}/include)
    set(ParMetis_DIR ${PARMETIS_INSTALL_PREFIX})
endif()

option(ENABLE_SUPERLU "Build the SuperLU solver package" TRUE)
if (ENABLE_SUPERLU)
  # SuperLU -- needs to go before HYPRE and PETSc because they needs it
  include(${SuperBuild_BUILD_FILES_DIR}/Build_SuperLU.cmake)
  append_set(Amanzi_TPL_CMAKE_ARGS
             -DSuperLU_DIR:FILEPATH=${TPL_INSTALL_PREFIX})

  # SuperLUDist
  include(${SuperBuild_BUILD_FILES_DIR}/Build_SuperLUDist.cmake)
  append_set(Amanzi_TPL_CMAKE_ARGS
    -DSuperLUDist_DIR:FILEPATH=${TPL_INSTALL_PREFIX})
endif()

# HYPRE
option(ENABLE_HYPRE "Build the HYPRE solver package" TRUE)
if (ENABLE_HYPRE)
  if (NOT ENABLE_XSDK)
    include(${SuperBuild_BUILD_FILES_DIR}/Build_HYPRE.cmake)
    set(HYPRE_INSTALL_PREFIX ${TPL_INSTALL_PREFIX})
    set(HYPRE_DIR ${HYPRE_INSTALL_PREFIX})
    set(HYPRE_INCLUDE_DIRS "${HYPRE_INSTALL_PREFIX}/include")
  else()
    find_spack_install_location(hypre HYPRE_INSTALL_PREFIX)
    message(STATUS "HYPRE_INSTALL_PREFIX = ${HYPRE_INSTALL_PREFIX}")
    set(HYPRE_DIR ${HYPRE_INSTALL_PREFIX})
    set(HYPRE_INCLUDE_DIRS "${HYPRE_INSTALL_PREFIX}/include")
    build_library_name(HYPRE HYPRE_LIBRARY APPEND_PATH ${HYPRE_INSTALL_PREFIX}/lib)
    set(HYPRE_LIBRARIES "${HYPRE_LIBRARY}")
  endif()
endif()

if (ENABLE_PETSC OR ENABLE_Structured)
  append_set(Amanzi_TPL_CMAKE_ARGS
             -DParMetis_DIR:FILEPATH=${ParMetis_DIR})
endif()

# PETSc
option(ENABLE_PETSC "Build the PETSc package" FALSE)
if (ENABLE_PETSC OR ENABLE_Structured)
  set(ENABLE_PETSC TRUE)
  if (NOT ENABLE_XSDK)
    include(${SuperBuild_BUILD_FILES_DIR}/Build_PETSc.cmake)
  else()
    find_spack_install_location(petsc PETSC_INSTALL_PREFIX)
    message(STATUS "PETSC_INSTALL_PREFIX = ${PETSC_INSTALL_PREFIX}")
    set(PETSC_DIR ${PETSC_INSTALL_PREFIX})
  endif()
endif()

# SEACAS
include(${SuperBuild_BUILD_FILES_DIR}/Build_SEACAS.cmake)

# Trilinos  
if (NOT ENABLE_XSDK)
  include(${SuperBuild_BUILD_FILES_DIR}/Build_Trilinos.cmake)
else()

  find_spack_install_location(trilinos Trilinos_INSTALL_PREFIX)
  message(STATUS "Trilinos_INSTALL_PREFIX = ${Trilinos_INSTALL_PREFIX}")

  find_spack_install_location(trilinos Zoltan_INSTALL_PREFIX)
  message(STATUS "Zoltan_INSTALL_PREFIX = ${Zoltan_INSTALL_PREFIX}")

endif()
append_set(Amanzi_TPL_CMAKE_ARGS
           -DENABLE_STK_Mesh:BOOL=${ENABLE_STK_Mesh}
           -DTrilinos_INSTALL_PREFIX:PATH=${Trilinos_INSTALL_PREFIX})

# MSTK
if (ENABLE_MSTK_Mesh)
  include(${SuperBuild_BUILD_FILES_DIR}/Build_MSTK.cmake) 
  append_set(Amanzi_TPL_CMAKE_ARGS
             -DENABLE_MSTK_Mesh:BOOL=${ENABLE_MSTK_Mesh}
             -DMSTK_INCLUDE_DIR:PATH=${MSTK_INCLUDE_DIR}
             -DMSTK_LIBRARY_DIR:PATH=${MSTK_LIBRARY_DIR})
endif()

# Silo
if (ENABLE_Silo)
  include(${SuperBuild_BUILD_FILES_DIR}/Build_Silo.cmake)
  append_set(Amanzi_TPL_CMAKE_ARGS
             -DENABLE_Silo:BOOL=${ENABLE_Silo}
             -DSilo_INCLUDE_DIR:PATH=${Silo_INCLUDE_DIR}
             -DSilo_LIBRARY_DIR:PATH=${Silo_LIBRARY_DIR})
endif()

# CLM
if (ENABLE_CLM)
  include(${SuperBuild_BUILD_FILES_DIR}/Build_CLM.cmake)
endif()

# Alquimia and chemical engines (PFloTran and CrunchTope)
if (NOT ENABLE_XSDK)
  
  option(ENABLE_ALQUIMIA "Build the Alquimia chemistry interface" FALSE)
  option(ENABLE_PFLOTRAN "Build the PFlotran chemistry package" FALSE)
  option(ENABLE_CRUNCHTOPE "Build the CrunchTope chemistry package" FALSE)

  # Alqumia depends on PFloTran
  if (ENABLE_ALQUIMIA)
    if (NOT ENABLE_PFLOTRAN)
      set(ENABLE_PFLOTRAN TRUE)
    endif()
  endif()

  if (ENABLE_PFLOTRAN)
    if (NOT ENABLE_PETSC)
      set(ENABLE_PETSC TRUE)
      include(${SuperBuild_BUILD_FILES_DIR}/Build_PETSc.cmake)
    endif()
    include(${SuperBuild_BUILD_FILES_DIR}/Build_PFloTran.cmake)
  endif()

  if (ENABLE_CRUNCHTOPE)
    if (NOT ENABLE_PETSC)
      set(ENABLE_PETSC TRUE)
      include(${SuperBuild_BUILD_FILES_DIR}/Build_PETSc.cmake)
    endif()
    include(${SuperBuild_BUILD_FILES_DIR}/Build_CrunchTope.cmake)
  endif()

  if (ENABLE_ALQUIMIA)
    include(${SuperBuild_BUILD_FILES_DIR}/Build_Alquimia.cmake)
  endif()

else()
  
  # pflotran
  find_spack_install_location(pflotran PFLOTRAN_INSTALL_PREFIX)
  message(STATUS "PFLOTRAN_INSTALL_PREFIX = ${PFLOTRAN_INSTALL_PREFIX}")
  build_library_name(pflotranchem PFLOTRAN_LIBRARIES APPEND_PATH ${PFLOTRAN_INSTALL_PREFIX}/lib)
  set(PFLOTRAN_INCLUDE_DIRS ${PFLOTRAN_INSTALL_PREFIX}/include)
  set(PFLOTRAN_DIR ${PFLOTRAN_INSTALL_PREFIX})
    
  # crunchtope
  if (ENABLE_CRUNCHTOPE)
    if (NOT ENABLE_PETSC)
      set(ENABLE_PETSC TRUE)
      #include(${SuperBuild_BUILD_FILES_DIR}/Build_PETSc.cmake)
    endif()
    include(${SuperBuild_BUILD_FILES_DIR}/Build_CrunchTope.cmake)
    build_library_name(crunchchem CRUNCHTOPE_LIBRARIES APPEND_PATH ${TPL_INSTALL_PREFIX}/lib)
    set(CRUNCHTOPE_INCLUDE_DIRS ${TPL_INSTALL_PREFIX}/include)
    set(CRUNCHTOPE_DIR ${TPL_INSTALL_PREFIX})
  endif()

  # alquimia
  find_spack_install_location(alquimia ALQUIMIA_INSTALL_PREFIX)
  message(STATUS "ALQUIMIA_INSTALL_PREFIX = ${ALQUIMIA_INSTALL_PREFIX}")
  build_library_name(alquimia_c ALQUIMIA_C_LIB APPEND_PATH ${ALQUIMIA_INSTALL_PREFIX}/lib)
  build_library_name(alquimia_cutils ALQUIMIA_CUTILS_LIB APPEND_PATH ${ALQUIMIA_INSTALL_PREFIX}/lib)
  build_library_name(alquimia_fortran ALQUIMIA_F_LIB APPEND_PATH ${ALQUIMIA_INSTALL_PREFIX}/lib)
  set(ALQUIMIA_INSTALL_PREFIX ${ALQUIMIA_INSTALL_PREFIX})

endif()

# Nanoflann is used by MSTK only
if (ENABLE_MSTK_Mesh)
  include(${SuperBuild_BUILD_FILES_DIR}/Build_Nanoflann.cmake)
  set(Nanoflann_INSTALL_PREFIX ${TPL_INSTALL_PREFIX})
  set(Nanoflann_INCLUDE_DIRS "${Nanoflann_INSTALL_PREFIX}/include")
endif()

# EXPRTK
include(${SuperBuild_BUILD_FILES_DIR}/Build_exprtk.cmake)


# ############################################################################ #
#  Final steps
# ############################################################################ #
# -- override compiler with nvccwrapper
if (ENABLE_KOKKOS_CUDA)
  set(NVCC_WRAPPER_DEFAULT_COMPILER ${CMAKE_CXX_COMPILER})
  set(CMAKE_CXX_COMPILER ${Trilinos_INSTALL_PREFIX}/bin/nvcc_wrapper)
else()
  set(NVCC_WRAPPER_DEFAULT_COMPILER "")
endif()

#  --- Create a CMake cache script
set(SuperBuild_Amanzi_CACHE_SCRIPT amanzi-tpl-config.cmake)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/${SuperBuild_Amanzi_CACHE_SCRIPT}.in
               ${SuperBuild_BINARY_DIR}/${SuperBuild_Amanzi_CACHE_SCRIPT}
               @ONLY)
install(FILES "${SuperBuild_BINARY_DIR}/${SuperBuild_Amanzi_CACHE_SCRIPT}"             
        DESTINATION ${TPL_INSTALL_PREFIX}/share/cmake)

# --- Create a Exodus split mesh script and add to the install target     
set(SuperBuild_EXOSPLIT_SCRIPT splitexo)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/${SuperBuild_EXOSPLIT_SCRIPT}.sh.in
               ${SuperBuild_BINARY_DIR}/${SuperBuild_EXOSPLIT_SCRIPT}
               @ONLY)
install(PROGRAMS "${SuperBuild_BINARY_DIR}/${SuperBuild_EXOSPLIT_SCRIPT}"             
         DESTINATION ${TPL_INSTALL_PREFIX}/bin)

#  --- Build the docs
option(BUILD_DOCS "Build the SuperBuild documentation" FALSE)
if (BUILD_DOCS)
  add_subdirectory(doc)
endif()


