#  -*- mode: cmake -*-

#
# Amanzi SuperBuild
#       Amanzi TPL Build Targets
cmake_minimum_required(VERSION 2.8.6)

if ("${CMAKE_PATCH_VERSION}" GREATER 3)
  set(ADJUST_POLICY "True")
endif()

if (${ADJUST_POLICY})
  cmake_policy(SET CMP0017 NEW)
endif()

# ############################################################################ #
# Define project name (SuperBuild_SOURCE_DIR and SuperBuild_BINARY_DIR)
# ############################################################################ #
project(SuperBuild)

# ############################################################################ #
# Define Amanzi directories 
# ############################################################################ #
file(TO_CMAKE_PATH ${SuperBuild_SOURCE_DIR}/../.. Amanzi_SOURCE_DIR)
file(TO_CMAKE_PATH ${SuperBuild_BINARY_DIR}/../.. Amanzi_BINARY_DIR)

# ############################################################################ #
# Update the CMake module path
# ############################################################################ #
set(Amanzi_MODULE_PATH
           ${Amanzi_SOURCE_DIR}/tools/cmake
           ${Amanzi_SOURCE_DIR}/tools/cmake/Utils
           ${Amanzi_SOURCE_DIR}/tools/cmake/Modules)

set(CMAKE_MODULE_PATH 
         ${CMAKE_MODULE_PATH} 
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
if ( NOT CMAKE_C_COMPILER )
  message(WARNING "C compiler not specified. CMake will guess!")
endif()
enable_language(C)

if ( NOT CMAKE_CXX_COMPILER )
  message(WARNING "C++ compiler not specified. CMake will guess!")
endif()
enable_language(CXX)

if ( NOT CMAKE_Fortran_COMPILER )
  message(WARNING "Fortran compiler not specified. CMake will guess!")
endif()
enable_language(Fortran)

include(FortranCInterface)
FortranCInterface_VERIFY()

include(DefineCompilerVersion)
define_compiler_version()

# ############################################################################ #
# Build options
# ############################################################################ #

message(STATUS "Setting build type")
set(SuperBuild_BUILD_TYPE_DFLT "Release")
if ( NOT CMAKE_BUILD_TYPE ) 
  set(CMAKE_BUILD_TYPE ${SuperBuild_BUILD_TYPE_DFLT})
endif()
message(STATUS "Setting build type -- ${CMAKE_BUILD_TYPE}")

# Download directory TPL_DOWNLOAD_DIR
set(TPL_DOWNLOAD_DIR_DFLT ${SuperBuild_BINARY_DIR}/Downloads)
if(NOT TPL_DOWNLOAD_DIR)
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
  check_download(TEST_FILE MD5SUMS TEST_URL http://software.lanl.gov/ascem/tpls)
endif()

# TPL install location
set(TPL_INSTALL_PREFIX_DFLT "${SuperBuild_BINARY_DIR}/external-projects")
if ( NOT TPL_INSTALL_PREFIX )
  set(TPL_INSTALL_PREFIX ${TPL_INSTALL_PREFIX_DFLT})
endif()
message(STATUS "Install TPLs in ${TPL_INSTALL_PREFIX}")

# Amanzi mesh build options
option(ENABLE_Structured   "Enable Amanzi structured mesh capability" TRUE)
option(ENABLE_Unstructured "Enable Amanzi unstructured mesh capability" TRUE)

option(ENABLE_STK_Mesh "Enable the Unstructured STK Mesh Toolkit" TRUE)
option(ENABLE_MSTK_Mesh "Enable the Unstructured MSTK Mesh Toolkit" TRUE)
set(ENABLE_MOAB_Mesh FALSE)

# OpenMP Search
option(ENABLE_OpenMP "Add OpenMP directives to the build" TRUE)
if ( ENABLE_OpenMP )
  find_package(OpenMP)
  if ( OPENMP_FOUND )
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  else()  
    message(SEND_ERROR "Failed to determine the OpenMP compiler flags. Will deactivate.")
    set(ENABLE_OpenMP FALSE CACHE BOOL "Enable OpenMP compile directives" FORCE)
  endif()
endif()


# Amanzi Test Suite
option(ENABLE_TESTS "Build the Amanzi Test Suite" TRUE)

# ############################################################################ #
# MPI Configuration
# ############################################################################ #

# MPI executable
if ( NOT MPI_EXEC )
  find_program(MPI_EXEC
               NAMES mpirun mpiexec aprun
               HINTS ENV MPI_ROOT ENV MPIROOT ENV MPI_PREFIX
               )
endif()

# Maximum number of MPI ranks
set(MPI_EXEC_MAX_NUMPROCS_DFLT 8)
if ( NOT MPI_EXEC_MAX_NUMPROCS )
  include(ProcessorCount)
  ProcessorCount(proc_count)
  if ( NOT proc_count EQUAL 0 )
    math(EXPR MPI_EXEC_MAX_NUMPROCS "${proc_count} * 2") 
    message(STATUS "Detected ${proc_count} processors and will set maximum to ${MPI_EXEC_MAX_NUMPROCS}")
  else()
    set(MPI_EXEC_MAX_NUMPROCS ${MPI_EXEC_MAX_NUMPROCS_DFLT})
  endif()
endif()  

# Number of MPI ranks flag
set(MPI_EXEC_NUMPROCS_FLAG_DFLT -n)
if(NOT MPI_EXEC_NUMPROCS_FLAG )
  set(MPI_EXEC_NUMPROCS_FLAG ${MPI_EXEC_NUMPROCS_FLAG_DFLT})
endif()

# ############################################################################ #
# Set common build compiler flags, build types and directories
# ############################################################################ #

# CMake compiler settings for any package built with CMake

# C language flags
set(Amanzi_CMAKE_C_COMPILER_ARGS
     -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
     -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
     -DCMAKE_C_FLAGS_DEBUG:STRING=${CMAKE_C_FLAGS_DEBUG}
     -DCMAKE_C_FLAGS_MINSIZEREL:STRING=${CMAKE_C_FLAGS_MINSIZEREL}
     -DCMAKE_C_FLAGS_RELEASE:STRING=${CMAKE_C_FLAGS_RELEASE}
     -DCMAKE_C_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_C_FLAGS_RELWITHDEBINFO})

# C++ language flags 
set(Amanzi_CMAKE_CXX_COMPILER_ARGS
     -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
     -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
     -DCMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}
     -DCMAKE_CXX_FLAGS_MINSIZEREL:STRING=${CMAKE_CXX_FLAGS_MINSIZEREL}
     -DCMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}
     -DCMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_CXX_FLAGS_RELWITHDEBINFO})

# Fortran language flags    
set(Amanzi_CMAKE_Fortran_COMPILER_ARGS
     -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
     -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS}
     -DCMAKE_Fortran_FLAGS_DEBUG:STRING=${CMAKE_Fortran_FLAGS_DEBUG}
     -DCMAKE_Fortran_FLAGS_MINSIZEREL:STRING=${CMAKE_Fortran_FLAGS_MINSIZEREL}
     -DCMAKE_Fortran_FLAGS_RELEASE:STRING=${CMAKE_Fortran_FLAGS_RELEASE}
     -DCMAKE_Fortran_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_Fortran_FLAGS_RELWITHDEBINFO})

# Link flags
set(Amanzi_CMAKE_LINKER_ARGS
     -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
     -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS})

# GNU configure static/shared
set(Amanzi_SHARED_SWITCH "--disable-shared")
if(BUILD_SHARED_LIBS)
  set(Amanzi_SHARED_SWITCH "--enable-shared")
endif()


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
##message(STATUS "\tC flags\t\t${Amanzi_COMMON_CFLAGS}")
#message(STATUS "\tC++ flags\t${Amanzi_COMMON_CXXFLAGS}")
#message(STATUS "\tFortran flags\t${Amanzi_COMMON_FCFLAGS}")

# ############################################################################ #
# Begin TPL builds 
# ############################################################################ #

# --- Initalize and set the TPL build files

# Create configure, build, install and test targets for each TPL
set_property(DIRECTORY PROPERTY
             EP_STEP_TARGETS download patch configure build install test)

# Include the TPL version information
include(${SuperBuild_SOURCE_DIR}/TPLVersions.cmake)

# Include the CMake module ExternalProject
include(ExternalProject)

# Include the macro that defines uniform build, install and Add_ExternalProject args
include(DefineExternalProjectArgs)

# Set of CMake configure arguments to pass to Amanzi once the TPLs are built
set(Amanzi_TPL_CMAKE_ARGS)

# --- Begin the build definitions

# ZLIB
include(${SuperBuild_BUILD_FILES_DIR}/Build_zlib.cmake) 

# CURL
include(${SuperBuild_BUILD_FILES_DIR}/Build_CURL.cmake)

# UnitTest
append_set(Amanzi_TPL_CMAKE_ARGS
           -DENABLE_TEST:BOOL=${ENABLE_TESTS})
if (ENABLE_TESTS)
  include(${SuperBuild_BUILD_FILES_DIR}/Build_UnitTest.cmake)
  append_set(Amanzi_TPL_CMAKE_ARGS 
             -DUnitTest_DIR:FILEPATH=${TPL_INSTALL_PREFIX})
endif()

# METIS  
include(${SuperBuild_BUILD_FILES_DIR}/Build_METIS.cmake) 
append_set(Amanzi_TPL_CMAKE_ARGS 
           -DMETIS_DIR:FILEPATH=${METIS_install_dir})
# ASCEM-IO  
#include(${SuperBuild_BUILD_FILES_DIR}/Build_ASCEMIO.cmake) 
#append_set(Amanzi_TPL_CMAKE_ARGS
#           -DASCEMIO_DIR:FILEPATH=${TPL_INSTALL_PREFIX})

# Now build the remaining TPLs

#DEBUG
#DEBUG# CCSE
#DEBUG# Temp option until the CMakeLists.txt file is fixed
#DEBUGset(CCSE_BL_SPACEDIM_DFLT 2)
#DEBUGif(NOT CCSE_BL_SPACEDIM )
#DEBUG  set(CCSE_BL_SPACEDIM ${CCSE_BL_SPACEDIM_DFLT})
#DEBUGendif()
#DEBUGinclude(${SuperBuild_BUILD_FILES_DIR}/Build_CCSE.cmake) 
#DEBUGappend_set(Amanzi_TPL_CMAKE_ARGS 
#DEBUG           -DCCSE_DIR:FILEPATH=${CCSE_install_dir}
#DEBUG           -DENABLE_MPI:BOOL=TRUE
#DEBUG           -DENABLE_OpenMP:BOOL=${ENABLE_OpenMP}
#DEBUG           -DAMANZI_PRECISION:STRING=DOUBLE
#DEBUG           -DAMANZI_SPACEDIM:INT=${CCSE_BL_SPACEDIM})
#DEBUG           
#DEBUG# CGNS
#DEBUGinclude(${SuperBuild_BUILD_FILES_DIR}/Build_CGNS.cmake) 
#DEBUGappend_set(Amanzi_TPL_CMAKE_ARGS 
#DEBUG          -DENABLE_CGNS:BOOL=TRUE
#DEBUG          -DCGNS_DIR:FILEPATH=${CGNS_install_dir})
#DEBUG
#DEBUG# Boost
#DEBUGinclude(${SuperBuild_BUILD_FILES_DIR}/Build_Boost.cmake) 
#DEBUG#TESTINGinclude(${SuperBuild_BUILD_FILES_DIR}/Build_BoostCmake.cmake) 
#DEBUGappend_set(Amanzi_TPL_CMAKE_ARGS
#DEBUG           -DBOOST_ROOT:FILEPATH=${Boost_install_dir}
#DEBUG	   -DBoost_USE_STATIC_LIBS:BOOL=TRUE
#DEBUG	   -DBoost_NO_SYSTEM_PATHS:BOOL=TRUE)
#DEBUG
#DEBUG# HDF5  
#DEBUGinclude(${SuperBuild_BUILD_FILES_DIR}/Build_HDF5.cmake) 
#DEBUG#BROKENinclude(${SuperBuild_BUILD_FILES_DIR}/Build_HDF5Cmake.cmake) 
#DEBUGappend_set(Amanzi_TPL_CMAKE_ARGS
#DEBUG           -DHDF5_DIR:FILEPATH=${HDF5_install_dir})
#DEBUG
#DEBUG
#DEBUG# NetCDF
#DEBUGinclude(${SuperBuild_BUILD_FILES_DIR}/Build_netCDF.cmake) 
#DEBUGappend_set(Amanzi_TPL_CMAKE_ARGS
#DEBUG           -DNetCDF_DIR:FILEPATH=${NETCDF_install_dir})
#DEBUG
#DEBUG# ExodusII
#DEBUGinclude(${SuperBuild_BUILD_FILES_DIR}/Build_ExodusII.cmake) 
#DEBUGappend_set(Amanzi_TPL_CMAKE_ARGS
#DEBUG           -DExodusII_DIR:FILEPATH=${ExodusII_install_dir})
#DEBUG
#DEBUG# MSTK
#DEBUGinclude(${SuperBuild_BUILD_FILES_DIR}/Build_MSTK.cmake) 
#DEBUGappend_set(Amanzi_TPL_CMAKE_ARGS
#DEBUG  -DENABLE_MSTK_Mesh:BOOL=${ENABLE_MSTK_Mesh}
#DEBUG          -DMSTK_INCLUDE_DIR:PATH=${MSTK_INCLUDE_DIR}
#DEBUG          -MSTK_LIBRARY_DIR:PATH=${MSTK_LIBRARY_DIR})
#DEBUG
#DEBUG# Trilinos  
#DEBUGinclude(${SuperBuild_BUILD_FILES_DIR}/Build_Trilinos.cmake)
#DEBUGappend_set(Amanzi_TPL_CMAKE_ARGS
#DEBUG           -DENABLE_STK_Mesh:BOOL=TRUE
#DEBUG           -DTrilinos_DIR:FILEPATH=${Trilinos_install_dir})
#DEBUG
#DEBUG# ############################################################################ #
#DEBUG# Create a CMake cache scipt
#DEBUG# ############################################################################ #
#DEBUGset(SuperBuild_Amanzi_CACHE_SCRIPT amanzi-tpl-config.cmake)
#DEBUGconfigure_file(${SuperBuild_BUILD_FILES_DIR}/${SuperBuild_Amanzi_CACHE_SCRIPT}.in
#DEBUG               ${SuperBuild_BINARY_DIR}/${SuperBuild_Amanzi_CACHE_SCRIPT}
#DEBUG               @ONLY)
#DEBUGinstall( FILES "${SuperBuild_BINARY_DIR}/${SuperBuild_Amanzi_CACHE_SCRIPT}"             
#DEBUG         DESTINATION ${TPL_INSTALL_PREFIX}/share/cmake
#DEBUG       )
#DEBUG
#DEBUG# ############################################################################ #
#DEBUG#  Final steps
#DEBUG# ############################################################################ #
#DEBUG
#DEBUG# Build the docs
#DEBUGadd_subdirectory(doc)
#DEBUG
#DEBUG# Build the Amanzi build script
#DEBUGoption(BUILD_Amanzi "Add Amanzi to the build stack" FALSE)
#DEBUGif ( BUILD_Amanzi )
#DEBUG
#DEBUG  set(Amanzi_CMAKE_ARGS ${Amanzi_CMAKE_COMPILER_ARGS})
#DEBUG
#DEBUG  # Enable the config report
#DEBUG  append_set(Amanzi_CMAKE_ARGS -DENABLE_Config_Report:BOOL=TRUE)
#DEBUG
#DEBUG  # Install path
#DEBUG  if (Amanzi_INSTALL_PREFIX)
#DEBUG    append_set(Amanzi_CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:FILEPATH=${Amanzi_INSTALL_PREFIX})
#DEBUG  endif()
#DEBUG
#DEBUG  # Build type
#DEBUG  if(CMAKE_BUILD_TYPE)
#DEBUG    append_set(Amanzi_CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})
#DEBUG  endif() 
#DEBUG
#DEBUG  # Disable shared libraries and executable unless user has requested shared libs
#DEBUG  set(Amanzi_SHARED_CMAKE_ARGS)
#DEBUG  if(BUILD_SHARED_LIBS)
#DEBUG    append_set(Amanzi_SHARED_CMAKE_ARGS -DBUILD_SHARED_LIBS:BOOL=TRUE)
#DEBUG  else()  
#DEBUG    append_set(Amanzi_SHARED_CMAKE_ARGS 
#DEBUG               -DBUILD_SHARED_LIBS:BOOL=FALSE
#DEBUG            -DPREFER_STATIC_LIBRARIES:BOOL=TRUE)
#DEBUG    #        -DBUILD_STATIC_EXECUTABLES:BOOL=TRUE)
#DEBUG  endif()      
#DEBUG
#DEBUG  append_set(Amanzi_CMAKE_ARGS ${Amanzi_SHARED_CMAKE_ARGS})
#DEBUG
#DEBUG  # MPI
#DEBUG  set(Amanzi_MPI_CMAKE_ARGS)
#DEBUG  if ( MPI_EXEC )
#DEBUG    append_set(Amanzi_MPI_CMAKE_ARGS -DMPI_EXEC:FILEPATH=${MPI_EXEC})
#DEBUG  else()  
#DEBUG    message(FATAL_ERROR "Enable Amanzi build requires MPI parallel run launch command\n"
#DEBUG                        "Please define the this binary with:\n"
#DEBUG                        "-DMPI_EXEC:FILEPATH=<binary>\n"
#DEBUG                        "and re-run cmake")
#DEBUG  endif()        
#DEBUG
#DEBUG  if ( MPI_EXEC_MAX_NUMPROCS )
#DEBUG    append_set(Amanzi_MPI_CMAKE_ARGS -DMPI_EXEC_MAX_NUMPROCS:STRING=${MPI_EXEC_NUMPROCS})
#DEBUG  else()  
#DEBUG    append_set(Amanzi_MPI_CMAKE_ARGS -DMPI_EXEC_MAX_NUMPROCS:STRING=4)
#DEBUG  endif()        
#DEBUG
#DEBUG  if ( MPI_EXEC_NUMPROCS_FLAG )
#DEBUG    append_set(Amanzi_MPI_CMAKE_ARGS -DMPI_EXEC_NUMPROCS_FLAG:STRING=${MPI_EXEC_NUMPROCS_FLAG})
#DEBUG  else()  
#DEBUG    append_set(Amanzi_MPI_CMAKE_ARGS -DMPI_EXEC_NUMPROCS_FLAG:STRING=-n)
#DEBUG  endif()        
#DEBUG
#DEBUG  append_set(Amanzi_CMAKE_ARGS ${Amanzi_MPI_CMAKE_ARGS})
#DEBUG  
#DEBUG  append_set(Amanzi_CMAKE_ARGS ${Amanzi_TPL_CMAKE_ARGS})
#DEBUG
#DEBUG
#DEBUG  append_set(Amanzi_CMAKE_ARGS
#DEBUG             -DENABLE_Unstructured:BOOL=${ENABLE_Unstructured}
#DEBUG             -DENABLE_Structured:BOOL=${ENABLE_Structured})
#DEBUG
#DEBUG  #print_variable(Amanzi_CMAKE_ARGS) 
#DEBUG
#DEBUG  if(NOT Amanzi_SOURCE_DIR)
#DEBUG    message(FATAL_ERROR "Enable Amanzi build requires an Amanzi source directory."
#DEBUG                        "Define this dorectory with:"
#DEBUG                        "-DAmanzi_SOURCE_DIR:FILEPAHT=<Amanzi source location>"
#DEBUG                        "and re-run cmake")
#DEBUG  endif()
#DEBUG  message(STATUS "Will build Amanzi source located in ${Amanzi_SOURCE_DIR}")
#DEBUG
#DEBUG  if (ENABLE_Structured AND NOT ENABLE_CCSE )
#DEBUG    message(FATAL_ERROR "The structured part of Amanzi requires CCSE. "
#DEBUG                        "Please enable CCSE with"
#DEBUG			"\n-D ENABLE_CCSE:BOOL=TRUE\n"
#DEBUG			"and rerun cmake."
#DEBUG			)
#DEBUG  endif()
#DEBUG
#DEBUG
#DEBUG  # Now add amanzi to the build stack
#DEBUG  include(ExternalProject)
#DEBUG
#DEBUG  set(Amanzi_DEPENDS exodusii hdf5 ascemio unittest boost mstk cgns trilinos)
#DEBUG  if ( ENABLE_Structured )
#DEBUG    append_set(Amanzi_DEPENDS ccse)
#DEBUG  endif()  
#DEBUG  ExternalProject_Add(amanzi
#DEBUG                      DEPENDS ${Amanzi_DEPENDS}
#DEBUG                      SOURCE_DIR ${Amanzi_SOURCE_DIR}
#DEBUG                      BINARY_DIR ${CMAKE_BINARY_DIR}/amanzi-build
#DEBUG                      DOWNLOAD_COMMAND ""
#DEBUG                      UPDATE_COMMAND ""
#DEBUG                      CMAKE_ARGS ${Amanzi_CMAKE_ARGS}
#DEBUG                      INSTALL_COMMAND "")
#DEBUGendif()        
#DEBUG               
