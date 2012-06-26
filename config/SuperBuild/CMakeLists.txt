#  -*- mode: cmake -*-

#
# Amanzi SuperBuild
#       Amanzi TPL Build Targets
cmake_minimum_required(VERSION 2.8.5)

if ("${CMAKE_PATCH_VERSION}" GREATER 3)
  set(ADJUST_POLICY "True")
endif()

if (${ADJUST_POLICY})
  cmake_policy(SET CMP0017 NEW)
endif()

# ############################################################################ #
# Compiler Definition Warnings
# ############################################################################ #
if ( NOT CMAKE_C_COMPILER )
  message(WARNING "C compiler not specified. CMake will guess!")
endif()

if ( NOT CMAKE_CXX_COMPILER )
  message(WARNING "C++ compiler not specified. CMake will guess!")
endif()

if ( NOT CMAKE_Fortran_COMPILER )
  message(WARNING "Fortran compiler not specified. CMake will guess!")
endif()

# ############################################################################ #
# Define project name (SuperBuild_SOURCE_DIR and SuperBuild_BINARY_DIR)
# ############################################################################ #
project(SuperBuild C CXX Fortran)

# ############################################################################ #
# Define Amanzi directories 
# ############################################################################ #
file(TO_CMAKE_PATH ${SuperBuild_SOURCE_DIR}/../.. Amanzi_SOURCE_DIR)

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
include(FortranCInterface)
FortranCInterface_VERIFY()

include(DefineCompilerVersion)
define_compiler_version()

include(CheckMPISourceCompiles)
check_mpi_source_compiles(MPI_WRAPPERS_IN_USE)

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
option(ENABLE_MOAB_Mesh "Enable the Unstructured MOAB Mesh Toolkit" FALSE)

# OpenMP Search
option(ENABLE_OpenMP "Add OpenMP directives to the build" TRUE)
if ( ENABLE_OpenMP )
  find_package(OpenMP)
  if ( OPENMP_FOUND )
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    if ( ${CMAKE_C_COMPILER_ID} STREQUAL "Intel" )
      set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -openmp") 
      set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -openmp")
    endif() 
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
               HINTS ENV MPI_ROOT ENV MPIROOT ENV MPI_PREFIX ENV MPI_HOME ENV MPIHOME
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
           -DMETIS_DIR:FILEPATH=${TPL_INSTALL_PREFIX})

# CCSE
set(CCSE_BL_SPACEDIM_DFLT 2)
if(NOT CCSE_BL_SPACEDIM )
  set(CCSE_BL_SPACEDIM ${CCSE_BL_SPACEDIM_DFLT})
endif()
include(${SuperBuild_BUILD_FILES_DIR}/Build_CCSE.cmake) 
append_set(Amanzi_TPL_CMAKE_ARGS 
           -DCCSE_DIR:FILEPATH=${TPL_INSTALL_PREFIX}
           -DENABLE_MPI:BOOL=TRUE
           -DENABLE_OpenMP:BOOL=${ENABLE_OpenMP}
           -DAMANZI_PRECISION:STRING=DOUBLE
           -DAMANZI_SPACEDIM:INT=${CCSE_BL_SPACEDIM})
           
# CGNS
include(${SuperBuild_BUILD_FILES_DIR}/Build_CGNS.cmake) 
append_set(Amanzi_TPL_CMAKE_ARGS 
          -DENABLE_CGNS:BOOL=TRUE
          -DCGNS_DIR:FILEPATH=${TPL_INSTALL_PREFIX})
# Boost
include(${SuperBuild_BUILD_FILES_DIR}/Build_Boost.cmake) 
#TESTINGinclude(${SuperBuild_BUILD_FILES_DIR}/Build_BoostCmake.cmake) 
append_set(Amanzi_TPL_CMAKE_ARGS
           -DBOOST_ROOT:FILEPATH=${TPL_INSTALL_PREFIX}
           -DBoost_USE_STATIC_LIBS:BOOL=TRUE
           -DBoost_NO_SYSTEM_PATHS:BOOL=TRUE)

# HDF5  
include(${SuperBuild_BUILD_FILES_DIR}/Build_HDF5.cmake) 
#BROKENinclude(${SuperBuild_BUILD_FILES_DIR}/Build_HDF5Cmake.cmake) 
append_set(Amanzi_TPL_CMAKE_ARGS
           -DHDF5_ROOT:FILEPATH=${TPL_INSTALL_PREFIX})

# ASCEM-IO  
include(${SuperBuild_BUILD_FILES_DIR}/Build_ASCEMIO.cmake) 
append_set(Amanzi_TPL_CMAKE_ARGS
           -DASCEMIO_DIR:FILEPATH=${TPL_INSTALL_PREFIX})

# NetCDF
include(${SuperBuild_BUILD_FILES_DIR}/Build_NetCDF.cmake) 
append_set(Amanzi_TPL_CMAKE_ARGS
           -DNetCDF_DIR:FILEPATH=${TPL_INSTALL_PREFIX})

# ExodusII
include(${SuperBuild_BUILD_FILES_DIR}/Build_ExodusII.cmake) 
append_set(Amanzi_TPL_CMAKE_ARGS
           -DExodusII_DIR:FILEPATH=${TPL_INSTALL_PREFIX})

# MSTK
if ( ENABLE_MSTK_Mesh )
  include(${SuperBuild_BUILD_FILES_DIR}/Build_MSTK.cmake) 
  append_set(Amanzi_TPL_CMAKE_ARGS
            -DENABLE_MSTK_Mesh:BOOL=${ENABLE_MSTK_Mesh}
            -DMSTK_INCLUDE_DIR:PATH=${MSTK_INCLUDE_DIR}
            -DMSTK_LIBRARY_DIR:PATH=${MSTK_LIBRARY_DIR})
endif()

# MOAB
if ( ENABLE_MOAB_Mesh )
  include(${SuperBuild_BUILD_FILES_DIR}/Build_MOAB.cmake)
  append_set(Amanzi_TPL_CMAKE_ARGS
             -DENABLE_MOAB_Mesh:BOOL=${ENABLE_MOAB_Mesh}
	     -DMOAB_DIR:PATH=${TPL_INSTALL_PREFIX})
endif()

# Trilinos  
include(${SuperBuild_BUILD_FILES_DIR}/Build_Trilinos.cmake)
append_set(Amanzi_TPL_CMAKE_ARGS
           -DENABLE_STK_Mesh:BOOL=${ENABLE_STK_Mesh}
           -DTrilinos_DIR:FILEPATH=${TPL_INSTALL_PREFIX})

# SEACAS
include(${SuperBuild_BUILD_FILES_DIR}/Build_SEACAS.cmake)


# ############################################################################ #
#  Final steps
# ############################################################################ #

#  --- Create a CMake cache scipt
set(SuperBuild_Amanzi_CACHE_SCRIPT amanzi-tpl-config.cmake)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/${SuperBuild_Amanzi_CACHE_SCRIPT}.in
               ${SuperBuild_BINARY_DIR}/${SuperBuild_Amanzi_CACHE_SCRIPT}
               @ONLY)
install( FILES "${SuperBuild_BINARY_DIR}/${SuperBuild_Amanzi_CACHE_SCRIPT}"             
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
if ( BUILD_DOCS )
  add_subdirectory(doc)
endif()

# Build the Amanzi build script
option(BUILD_Amanzi "Add Amanzi to the build stack" FALSE)
if ( BUILD_Amanzi )

  set(Amanzi_CMAKE_ARGS ${Amanzi_CMAKE_COMPILER_ARGS})

  # Enable the config report
  append_set(Amanzi_CMAKE_ARGS -DENABLE_Config_Report:BOOL=TRUE)

  # Install path
  if (Amanzi_INSTALL_PREFIX)
    append_set(Amanzi_CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:FILEPATH=${Amanzi_INSTALL_PREFIX})
  endif()

  # Build type
  if(CMAKE_BUILD_TYPE)
    append_set(Amanzi_CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})
  endif() 

  # Disable shared libraries and executable unless user has requested shared libs
  set(Amanzi_SHARED_CMAKE_ARGS)
  if(BUILD_SHARED_LIBS)
    append_set(Amanzi_SHARED_CMAKE_ARGS -DBUILD_SHARED_LIBS:BOOL=TRUE)
  else()  
    append_set(Amanzi_SHARED_CMAKE_ARGS 
               -DBUILD_SHARED_LIBS:BOOL=FALSE
            -DPREFER_STATIC_LIBRARIES:BOOL=TRUE)
    #        -DBUILD_STATIC_EXECUTABLES:BOOL=TRUE)
  endif()      

  append_set(Amanzi_CMAKE_ARGS ${Amanzi_SHARED_CMAKE_ARGS})

  # MPI
  set(Amanzi_MPI_CMAKE_ARGS)
  if ( MPI_EXEC )
    append_set(Amanzi_MPI_CMAKE_ARGS -DMPI_EXEC:FILEPATH=${MPI_EXEC})
  else()  
    message(FATAL_ERROR "Enable Amanzi build requires MPI parallel run launch command\n"
                        "Please define the this binary with:\n"
                        "-DMPI_EXEC:FILEPATH=<binary>\n"
                        "and re-run cmake")
  endif()        

  if ( MPI_EXEC_MAX_NUMPROCS )
    append_set(Amanzi_MPI_CMAKE_ARGS -DMPI_EXEC_MAX_NUMPROCS:STRING=${MPI_EXEC_NUMPROCS})
  else()  
    append_set(Amanzi_MPI_CMAKE_ARGS -DMPI_EXEC_MAX_NUMPROCS:STRING=4)
  endif()        

  if ( MPI_EXEC_NUMPROCS_FLAG )
    append_set(Amanzi_MPI_CMAKE_ARGS -DMPI_EXEC_NUMPROCS_FLAG:STRING=${MPI_EXEC_NUMPROCS_FLAG})
  else()  
    append_set(Amanzi_MPI_CMAKE_ARGS -DMPI_EXEC_NUMPROCS_FLAG:STRING=-n)
  endif()        

  append_set(Amanzi_CMAKE_ARGS ${Amanzi_MPI_CMAKE_ARGS})
  
  append_set(Amanzi_CMAKE_ARGS ${Amanzi_TPL_CMAKE_ARGS})


  append_set(Amanzi_CMAKE_ARGS
             -DENABLE_Unstructured:BOOL=${ENABLE_Unstructured}
             -DENABLE_Structured:BOOL=${ENABLE_Structured})

  #print_variable(Amanzi_CMAKE_ARGS) 

  if(NOT Amanzi_SOURCE_DIR)
    message(FATAL_ERROR "Enable Amanzi build requires an Amanzi source directory."
                        "Define this dorectory with:"
                        "-DAmanzi_SOURCE_DIR:FILEPATH=<Amanzi source location>"
                        "and re-run cmake")
  endif()
  message(STATUS "Will build Amanzi source located in ${Amanzi_SOURCE_DIR}")


  # Now add amanzi to the build stack
  include(ExternalProject)

  set(Amanzi_DEPENDS 
                  ${ExodusII_BUILD_TARGET}
		  ${HDF5_BUILD_TARGET}
		  ${ASCEMIO_BUILD_TARGET}
		  ${UnitTest_BUILD_TARGET}
		  ${Boost_BUILD_TARGET}
		  ${MSTK_BUILD_TARGET}
		  ${MOAB_BUILD_TARGET}
		  ${CGNS_BUILD_TARGET}
		  ${CCSE_BUILD_TARGET}
		  ${Trilinos_BUILD_TARGET})

  if ( ENABLE_Structured )
    list(APPEND Amanzi_DEPENDS ${CCSE_BUILD_TARGET})
  endif()  
  ExternalProject_Add(amanzi
                      DEPENDS ${Amanzi_DEPENDS}
		      PREFIX ${CMAKE_BINARY_DIR}/Amanzi
                      SOURCE_DIR ${Amanzi_SOURCE_DIR}
                      DOWNLOAD_COMMAND ""
                      UPDATE_COMMAND ""
                      CMAKE_ARGS
		                ${Amanzi_CMAKE_ARGS}
				${Amanzi_CMAKE_C_COMPILER_ARGS}
				${Amanzi_CMAKE_CXX_COMPILER_ARGS}
				${Amanzi_CMAKE_Fortran_COMPILER_ARGS}
                      INSTALL_DIR ${Amanzi_INSTALL_PREFIX} )
endif()        
               
