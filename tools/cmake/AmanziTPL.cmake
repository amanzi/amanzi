# -*- mode: cmake -*-
# 
# Amanzi Third Party Library (TPL) Definitions
#

# Standard CMake modules see CMAKE_ROOT/Modules
include(FeatureSummary)

message(STATUS "")
message(STATUS ">>>>>>>> AmanziTPL.cmake")

set(AMANZI_TPLS_VERSION_CMAKE "${CMAKE_SOURCE_DIR}/config/SuperBuild/TPLVersions.cmake")
file(STRINGS ${AMANZI_TPLS_VERSION_CMAKE} AMANZI_TPLS_VERSION_CMAKE_lines)

FOREACH(line ${AMANZI_TPLS_VERSION_CMAKE_lines})
  IF ( ${line} MATCHES "AMANZI_TPLS_VERSION " )
    # message(STATUS "\t >>>>>> Skipping ....")
  ELSEIF ( ${line} MATCHES "AMANZI_TPLS_VERSION_MAJOR" )
    STRING(REGEX REPLACE ".*AMANZI_TPLS_VERSION_MAJOR ([0-9]+).*" "\\1" AMANZI_TPLS_VERSION_MAJOR_SRC "${line}")
    # message(STATUS "\t >>>>>> ${AMANZI_TPLS_VERSION_MAJOR_SRC}")
  ELSEIF ( ${line} MATCHES "AMANZI_TPLS_VERSION_MINOR" )
    STRING(REGEX REPLACE ".*AMANZI_TPLS_VERSION_MINOR ([0-9]+).*" "\\1" AMANZI_TPLS_VERSION_MINOR_SRC "${line}")
    # message(STATUS "\t >>>>>> ${AMANZI_TPLS_VERSION_MINOR_SRC}")
  ELSEIF ( ${line} MATCHES "AMANZI_TPLS_VERSION_PATCH" )
    STRING(REGEX REPLACE ".*AMANZI_TPLS_VERSION_PATCH ([0-9]+).*" "\\1" AMANZI_TPLS_VERSION_PATCH_SRC "${line}")
    # message(STATUS "\t >>>>>> ${AMANZI_TPLS_VERSION_PATCH_SRC}")
  ENDIF()
ENDFOREACH()

set(AMANZI_TPLS_VERSION_REQUIRED ${AMANZI_TPLS_VERSION_MAJOR_SRC}.${AMANZI_TPLS_VERSION_MINOR_SRC}.${AMANZI_TPLS_VERSION_PATCH_SRC})

# Verify TPL compatibility
# - should be able to report the path to the TPLs install, need to write it in cache.
message(STATUS "\t >>>>>> Amanzi Source Path:       ${CMAKE_SOURCE_DIR}")
message(STATUS "\t >>>>>  Amanzi Required Version:  ${AMANZI_TPLS_VERSION_REQUIRED}")
message(STATUS "\t >>>>>  Amanzi TPL Version:       ${AMANZI_TPLS_VERSION}\n")


if (NOT ${AMANZI_TPLS_VERSION} STREQUAL ${AMANZI_TPLS_VERSION_REQUIRED}) 
  message(WARNING "TPL version does not match the required version.")
endif()

# List of directories for step 1 in CMAKe's search path for command find_library(). 
# By default it is empty. We set it to non-standard location of MPI. 
# Probably, it would be appropriate to add also the path to TPL libraties
set(CMAKE_LIBRARY_PATH ${MPI_PREFIX}/lib)


# Amanzi CMake modules see <root source>/tools/cmake
include(CheckMPISourceCompiles)
include(TrilinosMacros)
include(PrintVariable)
include(AddImportedLibrary)

message(STATUS "CMAKE information")
message(STATUS "\tCMAKE_SYSTEM_LIBRARY_PATH = ${CMAKE_SYSTEM_LIBRARY_PATH}")
message(STATUS "\tCMAKE_LIBRARY_PATH = ${CMAKE_LIBRARY_PATH}")


##############################################################################
# ------------------------ Required Libraries -------------------------------#
##############################################################################

##############################################################################
# MPI
##############################################################################
check_mpi_source_compiles(MPI_WRAPPERS_IN_USE)

if ( NOT MPI_WRAPPERS_IN_USE )
  message(WARNING "At this time, Amanzi must be compiled with MPI wrappers."
                  " Build will likely fail. Please define CMAKE_*_COMPILER"
		  " parameters as MPI compiler wrappers and re-run cmake.")
endif()

##############################################################################
# ZLIB
##############################################################################

set (ZLIB_ROOT ${ZLIB_DIR})
find_package(ZLIB REQUIRED)

message(STATUS "Zlib Package information")
message(STATUS "\tZLIB_DIR          = ${ZLIB_DIR}")
message(STATUS "\tZLIB_INCLUDE_DIR  = ${ZLIB_INCLUDE_DIR}")
message(STATUS "\tZLIB_INCLUDE_DIRS = ${ZLIB_INCLUDE_DIRS}")
message(STATUS "\tZLIB_LIBRARY_DIR  = ${ZLIB_LIBRARY_DIR}")
message(STATUS "\tZLIB_LIBRARY      = ${ZLIB_LIBRARY}")
message(STATUS "\tZLIB_LIBRARIES    = ${ZLIB_LIBRARIES}")
message(STATUS "")

##############################################################################
# Boost
##############################################################################

# CMake 2.8.6 FindBoost stops at version 1.46
# Add more versions to the search see cmake --help-module FindBoost for
# more information.
set(Boost_ADDITIONAL_VERSIONS 
    1.47 1.47.0
    1.48 1.48.0
    1.49 1.49.0
    1.50 1.50.0
    1.51 1.51.0
    1.52 1.52.0
    1.53 1.53.0
    1.54 1.55.0)
find_package(Boost COMPONENTS system filesystem program_options regex REQUIRED)
set_package_properties(Boost PROPERTIES
                 DESCRIPTION "C++ Extension library"
                 URL "http://www.boost.org"
                 PURPOSE "Required by the MPC")

if (Boost_VERSION) 
  if ( ${Boost_VERSION} VERSION_LESS 1.46 )
    message(WARNING "Found Boost version ${Boost_VERSION} which"
                    " is older than the supported (1.46) version.")
  endif()

  # The Boost filesystem library changed and deprecated some functions.
  # This define should be used when packages include boost/filesystem.hpp
  # and packages any of these new or deprecated functions.
  # The change from version 2 to 3 occurred with the 1.49 Boost release.
  # Please refer to the online documentation at www.boost.org.
  if ( "${Boost_VERSION}" VERSION_LESS "1.34" )
    set(Boost_FILESYSTEM_DEFINES "BOOST_FILESYSTEM_VERSION=1")
  elseif ( "${Boost_VERSION}" VERSION_LESS "1.49" )
    set(Boost_FILESYSTEM_DEFINES "BOOST_FILESYSTEM_VERSION=2")
  else()
    set(Boost_FILESYSTEM_DEFINES "BOOST_FILESYSTEM_VERSION=3")
  endif()  
endif()


##############################################################################
# HDF5 - http://www.hdfgroup.org/HDF5/
##############################################################################
if (BUILD_SHARED_LIBS)
  set(HDF5_USE_STATIC_LIBRARIES OFF)
else()    
  set(HDF5_USE_STATIC_LIBRARIES ON)
endif()
  
find_package(HDF5 1.10.0 REQUIRED COMPONENTS C HL)

if (NOT HDF5_IS_PARALLEL) 
    message(WARNING "The HDF5 installation found in ${HDF5_DIR} is not "
                    "a parallel build. At this time, this installation "
                    "is compatible with other TPLs. Soon Amanzi will "
                    "require a parallel enabled HDF5. Please update your "
                    "HDF5 installation to include MPI I/O symbols"
           )            
endif(NOT HDF5_IS_PARALLEL)
set_package_properties(HDF5 PROPERTIES
                       DESCRIPTION "I/O library that creates HDF5 formatted files"
                       URL "http://www.hdfgroup.org/HDF5"
                       PURPOSE "Required library for several components in Amanzi"
                      )

set(HDF5_LIBRARIES ${HDF5_HL_LIBRARIES} ${HDF5_C_LIBRARIES} ${ZLIB_LIBRARIES} m dl)

message(STATUS "HDF5 Package information")
message(STATUS "\tHDF5_INCLUDE_DIR  = ${HDF5_INCLUDE_DIR}")
message(STATUS "\tHDF5_INCLUDE_DIRS = ${HDF5_INCLUDE_DIRS}")
message(STATUS "\tHDF5_LIBRARY_DIR  = ${HDF5_LIBRARY_DIR}")
message(STATUS "\tHDF5_LIBRARY      = ${HDF5_LIBRARY}")
message(STATUS "\tHDF5_LIBRARIES    = ${HDF5_LIBRARIES}")
message(STATUS "\tHDF5_HL_LIBRARIES = ${HDF5_HL_LIBRARIES}")
message(STATUS "")

##############################################################################
# Trilinos http://trilinos.sandia.gov
##############################################################################

# This command alters Trilinos_DIR. If it finds the configuration file
# Trilinos_DIR is set to the path the configuration file was found.
if ( NOT Trilinos_INSTALL_PREFIX )
  message(WARNING "Use Trilinos_INSTALL_PREFIX"
                  " to define the Trilinos installation location"
		  "\n-DTrilinos_INSTALL_PREFIX:PATH=<trilinos directory>\n")
endif()
set(Trilinos_MINIMUM_VERSION 13.0.0)
message(STATUS "Searching for Trilinos at: ${Trilinos_INSTALL_PREFIX}")
find_package(Trilinos ${Trilinos_MINIMUM_VERSION} REQUIRED
             PATHS ${Trilinos_INSTALL_PREFIX}
             PATH_SUFFIXES include)
            
if (Trilinos_FOUND)
  message(STATUS "Trilinos Package information")
  message(STATUS "\tTrilinos_VERSION      = ${Trilinos_VERSION}")
  message(STATUS "\tTrilinos_DIR          = ${Trilinos_DIR}")
  message(STATUS "\tTrilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
  message(STATUS "")

  trilinos_package_enabled_tpls(Trilinos)           

  if ("${Trilinos_VERSION}" VERSION_LESS ${Trilinos_MINIMUM_VERSION}) 
    message(FATAL_ERROR "Trilinos version ${Trilinos_VERSION} is not sufficient."
                        " Amanzi requires at least version ${Trilinos_MINIMUM_VERSION}")
  endif()

  # Amanzi uses Epetra and Teuchos utils throughout the code. 
  # This find_package call defines Epetra_* variables.
  # Amanzi developers should use these variables
  # for libraries that ONLY use Epetra/Teuchos and avoid
  # using the ALL POWERFUL(TM) Trilinos_LIBRARIES.
  # When/If we create wrappers, using this variable
  # will make that transition easier.

  # Verify that the Trilinos found has all the required packages
  # Amanzi needs. This list changes for structured and unstructured
  # mesh capabilities.
  # List of required Trilinos packages
  # Teuchos - general purpose toolkit, used through the code
  # Epetra  - distributed data objects
  # NOX     - nonlinear solver (Unstructured ONLY)
  # ML      - multilevel preconditioner (Unstructured ONLY)
  set(Trilinos_REQUIRED_PACKAGE_LIST Teuchos Epetra) 
  if (ENABLE_Unstructured)
    # list(APPEND Trilinos_REQUIRED_PACKAGE_LIST NOX ML Amesos2 MueLu)
    list(APPEND Trilinos_REQUIRED_PACKAGE_LIST NOX ML Amesos2)
  endif()

  foreach(tri_package ${Trilinos_REQUIRED_PACKAGE_LIST})
    find_package(${tri_package} REQUIRED
                 NO_MODULE
                 HINTS ${Trilinos_INSTALL_PREFIX}
                 PATH_SUFFIXES include lib)
    trilinos_package_enabled_tpls(${tri_package})
    message(STATUS "\t${tri_package}_DIR       = ${${tri_package}_DIR}")
    message(STATUS "\t${tri_package}_LIBRARIES = ${${tri_package}_LIBRARIES}")
    string(TOLOWER ${tri_package} _tri_package)
    print_link_libraries(${_tri_package})
    message(STATUS "")

    # Update the <PACKAGE>_INCLUDE_DIRS variable 
    foreach( _inc ${${tri_package}_TPL_INCLUDE_DIRS})
      list(APPEND ${tri_package}_INCLUDE_DIRS "${_inc}")
    endforeach()
  endforeach()

  # Zoltan - required by MSTK mesh class 
  if (ENABLE_MESH_MSTK)
    find_package(Zoltan
                 NO_MODULE
                 HINTS ${Trilinos_INSTALL_PREFIX}
                 PATH_SUFFIXES include lib)
    if (Zoltan_FOUND)
      message(STATUS "\tZoltan_DIR          = ${Zoltan_DIR}")
      message(STATUS "\tZoltan_INCLUDE_DIRS = ${Zoltan_INCLUDE_DIRS}")
      message(STATUS "\tZoltan_LIBRARIES    = ${Zoltan_LIBRARIES}")
      message(STATUS "")
      trilinos_package_enabled_tpls(Zoltan)
      foreach( _inc "${ZOLTAN_TPL_INCLUDE_DIRS}")
        list(APPEND ZOLTAN_INCLUDE_DIRS "${_inc}")
      endforeach()
    else()  
      message(WARNING "Could not locate Zoltan in ${Trilinos_DIR}. Will not enable MESH_MSTK")
      set(ENABLE_MESH_MSTK OFF CACHE BOOL "Disable MSTK Mesh capability" FORCE)
    endif()
  endif()

  option(ENABLE_HYPRE "Enable HYPRE APIs in flow" ON)
  if (ENABLE_HYPRE)
    # Ifpack - preconditioner package that serves as a wrapper for HYPRE
    find_package(Ifpack 
                 NO_MODULE
                 HINTS ${Trilinos_INSTALL_PREFIX}
                 PATH_SUFFIXES include lib
                 )
    if (Ifpack_FOUND)
      message(STATUS "\tIfpack_DIR          = ${Ifpack_DIR}")
      message(STATUS "\tIfpack_DIR          = ${Ifpack_DIR}")
      message(STATUS "\tIfpack_INCLUDE_DIRS = ${Ifpack_INCLUDE_DIRS}")
      message(STATUS "\tIfpack_LIBRARIES    = ${Ifpack_LIBRARIES}")
      print_link_libraries(ifpack)
      message(STATUS "")

      trilinos_package_enabled_tpls(Ifpack)

      foreach( _inc "${Ifpack_TPL_INCLUDE_DIRS}")
        list(APPEND Ifpack_INCLUDE_DIRS "${_inc}")
      endforeach()
    else()
      message(SEND_ERROR "Trilinos in ${Trilinos_DIR} does not have the Ifpack package")
    endif()

    if (NOT Ifpack_ENABLE_HYPRE AND NOT ENABLE_KOKKOS)
      message(WARNING "ENABLE_HYPRE requires the Trilinos package Ifpack with enabled HYPRE."
                      " Deactivating HYPRE APIs")
      set(ENABLE_HYPRE OFF CACHE BOOL "Disable the HYPRE APIs" FORCE)                 
    endif()
  endif()
       
  # Now update the Trilinos_LIBRARIES and INCLUDE_DIRS
  foreach( _inc "${Trilinos_TPL_INCLUDE_DIRS}")
    list(APPEND Trilinos_INCLUDE_DIRS "${_inc}")
    list(REMOVE_DUPLICATES Trilinos_INCLUDE_DIRS)
  endforeach()
else()
  message(FATAL_ERROR "Can not locate Trilinos configuration file\n"
                      " Please define the location of your Trilinos installation\n"
                      "using -D Trilinos_DIR:FILEPATH=<install path>\n")
endif()    


##############################################################################
# HYPRE and its dependencies
##############################################################################

if (ENABLE_HYPRE)

  find_package(HYPRE)

  if (HYPRE_FOUND)
    message(STATUS "HYPRE Package information")
    message(STATUS "\tHYPRE_VERSION      = ${HYPRE_VERSION}")
    message(STATUS "\tHYPRE_INCLUDE_DIRS = ${HYPRE_INCLUDE_DIRS}")
    message(STATUS "\tHYPRE_LIBRARY_DIR  = ${HYPRE_LIBRARY_DIR}")
    message(STATUS "\tHYPRE_LIBRARY      = ${HYPRE_LIBRARY}")
    message(STATUS "\tHYPRE_LIBRARIES    = ${HYPRE_LIBRARIES}")
    print_link_libraries(${HYPRE_LIBRARY})
    message(STATUS "")
  else()
    message(FATAL_ERROR "Can not locate HYPRE library and/or include\n")
  endif()
endif()


##############################################################################
# NetCDF - http://www.unidata.ucar.edu/software/netcdf/
##############################################################################

find_package(NetCDF REQUIRED)
set_package_properties(NetCDF
                 PROPERTIES
                 DESCRIPTION "Network Common Data Format (NetCDF)"
                 URL "http://www.unidata.ucar.edu/software/netcdf/"
                 PURPOSE "Required by ExodusII library")

message(STATUS "NetCDF Package information")
message(STATUS "\tNetCDF_VERSION      = ${NetCDF_VERSION}")
message(STATUS "\tNetCDF_INCLUDE_DIR  = ${NetCDF_INCLUDE_DIR}")
message(STATUS "\tNetCDF_INCLUDE_DIRS = ${NetCDF_INCLUDE_DIRS}")
message(STATUS "\tNetCDF_LIBRARY_DIR  = ${NetCDF_LIBRARY_DIR}")
message(STATUS "\tNetCDF_C_LIBRARY    = ${NetCDF_C_LIBRARY}")
message(STATUS "\tNetCDF_C_LIBRARIES  = ${NetCDF_C_LIBRARIES}")
message(STATUS "\tNetCDF_CXX_LIBRARIES  = ${NetCDF_CXX_LIBRARIES}")
print_link_libraries(${NetCDF_C_LIBRARY})
message(STATUS "")


##############################################################################
# XERCES-C - http://http://xerces.apache.org/xerces-c/
##############################################################################
find_package(XERCES REQUIRED)
set_package_properties(XERCES
                 PROPERTIES
                 PURPOSE "Validating XML parser")

if (XERCES_FOUND)
  message(STATUS "XERCES Package information")
  message(STATUS "\tXERCES_INCLUDE_DIR  = ${XERCES_INCLUDE_DIR}")
  message(STATUS "\tXERCES_INCLUDE_DIRS = ${XERCES_INCLUDE_DIRS}")
  message(STATUS "\tXERCES_LIBRARY_DIR  = ${XERCES_LIBRARY_DIR}")
  message(STATUS "\tXERCES_LIBRARY      = ${XERCES_LIBRARY}")
  message(STATUS "\tXERCES_LIBRARIES    = ${XERCES_LIBRARIES}")
  message(STATUS "\tXERCES_ICU_LIBRARIES= ${XERCES_ICU_LIBRARIES}")
  message(STATUS "")
endif()


##############################################################################
# SEACAS
##############################################################################
find_package(SEACAS REQUIRED)
set_package_properties(SEACAS
                 PROPERTIES
                 DESCRIPTION "FE analysis tools with ExodusII mesh (SEACAS)"
                 URL "https://github.com/gsjaardema/seacas"
                 PURPOSE "Required by mesh library")

set(SEACAS_LIBRARIES ${SEACAS_LIBRARIES} ${NetCDF_C_LIBRARIES} ${HDF5_LIBRARIES})

message(STATUS "SEACAS Package information")
message(STATUS "\tSEACAS_INCLUDE_DIR  = ${SEACAS_INCLUDE_DIR}")
message(STATUS "\tSEACAS_INCLUDE_DIRS = ${SEACAS_INCLUDE_DIRS}")
message(STATUS "\tSEACAS_LIBRARY_DIR  = ${SEACAS_LIBRARY_DIR}")
message(STATUS "\tSEACAS_LIBRARY      = ${SEACAS_LIBRARY}")
message(STATUS "\tSEACAS_LIBRARIES    = ${SEACAS_LIBRARIES}")
print_link_libraries(${SEACAS_LIBRARY})
message(STATUS "")


##############################################################################
# CCSE - http://ccse.lbl.gov/Software/ccse_core.html
##############################################################################
if (ENABLE_Structured)
  find_package(CCSE REQUIRED)
  set_package_properties(CCSE
                   PROPERTIES
                   PURPOSE "CCSE BoxLib softare library required for structured grid")

  if (CCSE_FOUND)
    message(STATUS "CCSE Package information")
    message(STATUS "\tCCSE_INCLUDE_DIR  = ${CCSE_INCLUDE_DIR}")
    message(STATUS "\tCCSE_INCLUDE_DIRS = ${CCSE_INCLUDE_DIRS}")
    message(STATUS "\tCCSE_LIBRARY_DIR  = ${CCSE_LIBRARY_DIR}")
    message(STATUS "\tCCSE_LIBRARY      = ${CCSE_LIBRARY}")
    message(STATUS "\tCCSE_LIBRARIES    = ${CCSE_LIBRARIES}")
    print_link_libraries(${CCSE_LIBRARY})
    message(STATUS "")
  endif()
endif()


##############################################################################
############################ Option Processing ###############################
##############################################################################


##############################################################################
#---------------------------- Mesh Frameworks -------------------------------#
##############################################################################

# Enable ALL possible mesh frameworks
#option(ENABLE_MESH_ALL "Build all Amanzi mesh frameworks" OFF)
#if(ENABLE_MESH_ALL)
#   set(ENABLE_MESH_MOAB ON)
#   set(ENABLE_MESH_MSTK ON)
#endif()    
#set_feature_info(MESH_ALL
#                 ENABLE_MESH_ALL
#                 "Build all available mesh frameworks"
#                )    


##############################################################################
# MOAB - svn co https://svn.mcs.anl.gov/repos/ITAPS/MOAB/trunk MOAB
##############################################################################
option(ENABLE_MESH_MOAB "Build Amanzi with the MOAB mesh framework" OFF)
add_feature_info(MESH_MOAB
                 ENABLE_MESH_MOAB
                 "A Mesh-Oriented datABase"
                 )
if (ENABLE_MESH_MOAB)
  find_package(MOAB REQUIRED)

  if (MOAB_FOUND)
    message(STATUS "MOAB Package information")
    message(STATUS "\tMOAB_VERSION      = ${MOAB_VERSION}")
    message(STATUS "\tMOAB_INCLUDE_DIR  = ${MOAB_INCLUDE_DIR}")
    message(STATUS "\tMOAB_INCLUDE_DIRS = ${MOAB_INCLUDE_DIRS}")
    message(STATUS "\tMOAB_LIBRARY_DIR  = ${MOAB_LIBRARY_DIR}")
    message(STATUS "\tMOAB_LIBRARY      = ${MOAB_LIBRARY}")
    message(STATUS "\tMOAB_LIBRARIES    = ${MOAB_LIBRARIES}")
    print_link_libraries(${MOAB_LIBRARY})
    message(STATUS "")
  endif() 
endif()


##############################################################################
# MSTK - https://software.lanl.gov/MeshTools/trac/raw-attachment/wiki/WikiStart/mstk-1.80.tar.gz
##############################################################################
option(ENABLE_MESH_MSTK "Build Amanzi with the MSTK mesh framework" OFF)
add_feature_info(MESH_MSTK
                 ENABLE_MESH_MSTK
                 "A mesh framework"
                 )
if (ENABLE_MESH_MSTK)
  find_package(MSTK REQUIRED)

  if (MSTK_FOUND)
    message(STATUS "MSTK Package information")
    message(STATUS "\tMSTK_VERSION      = ${MSTK_VERSION}")
    message(STATUS "\tMSTK_INCLUDE_DIR  = ${MSTK_INCLUDE_DIR}")
    message(STATUS "\tMSTK_INCLUDE_DIRS = ${MSTK_INCLUDE_DIRS}")
    message(STATUS "\tMSTK_LIBRARY_DIR  = ${MSTK_LIBRARY_DIR}")
    message(STATUS "\tMSTK_LIBRARY      = ${MSTK_LIBRARY}")
    message(STATUS "\tMSTK_LIBRARIES    = ${MSTK_LIBRARIES}")
    print_link_libraries(${MSTK_LIBRARY})
    message(STATUS "")
  endif() 
endif() 


##############################################################################
# Silo
##############################################################################
option(ENABLE_Silo "Build Amanzi with Silo output options" OFF)
if (ENABLE_Silo)
  find_package(Silo REQUIRED)

  if (Silo_FOUND)
    message(STATUS "Silo Package information")
    message(STATUS "\tSilo_INCLUDE_DIR  = ${Silo_INCLUDE_DIR}")
    message(STATUS "\tSilo_INCLUDE_DIRS = ${Silo_INCLUDE_DIRS}")
    message(STATUS "\tSilo_LIBRARY_DIR  = ${Silo_LIBRARY_DIR}")
    message(STATUS "\tSilo_LIBRARY      = ${Silo_LIBRARY}")
    message(STATUS "\tSilo_LIBRARIES    = ${Silo_LIBRARIES}")
    message(STATUS "")
  endif() 
endif() 


##############################################################################
#-------------------------- Optional Libraries ------------------------------#
##############################################################################


##############################################################################
# ASCEMIO - http://www.cgns.sourceforge.net/
##############################################################################
option(ENABLE_ASCEMIO  "Build Amanzi output library with ASCEM-IO parallelIO" OFF)
set_package_properties(ASCEMIO
                       PROPERTIES
                       DESCRIPTION "ASCEM-IO Scalable Parallel I/O module for Environmental Management Applications"
                       URL "http://ascem-io.secure-water.org"
                       PURPOSE "Required to produce VisIt files in parallel"
                      )
#if (ENABLE_ASCEMIO)
if (ENABLE_Unstructured)
  find_package(ASCEMIO REQUIRED)
else()
  find_package(ASCEMIO)
endif() 

if (ASCEMIO_FOUND)
  message(STATUS "ASCEMIO Package information")
  message(STATUS "\tASCEMIO_INCLUDE_DIR  = ${ASCEMIO_INCLUDE_DIR}")
  message(STATUS "\tASCEMIO_INCLUDE_DIRS = ${ASCEMIO_INCLUDE_DIRS}")
  message(STATUS "\tASCEMIO_LIBRARY_DIR  = ${ASCEMIO_LIBRARY_DIR}")
  message(STATUS "\tASCEMIO_LIBRARY      = ${ASCEMIO_LIBRARY}")
  message(STATUS "\tASCEMIO_LIBRARIES    = ${ASCEMIO_LIBRARIES}")
  message(STATUS "")
endif() 


##############################################################################
# UnitTest++ - http://unittest-cpp.sourceforge.net/
##############################################################################
option(ENABLE_UnitTest "Build Amanzi unit tests. Requires UnitTest++" ON)
set_package_properties(UnitTest
                       PROPERTIES
                       DESCRIPTION "C++ unit test framework"
                       URL "http://unittest-cpp.sourceforge.net/"
                      )
if (ENABLE_UnitTest)
  find_package(UnitTest)

  if (UNITTEST_FOUND)
    message(STATUS "UnitTest Package information")
    message(STATUS "\tUnitTest_INCLUDE_DIR  = ${UnitTest_INCLUDE_DIR}")
    message(STATUS "\tUnitTest_INCLUDE_DIRS = ${UnitTest_INCLUDE_DIRS}")
    message(STATUS "\tUnitTest_LIBRARY_DIR  = ${UnitTest_LIBRARY_DIR}")
    message(STATUS "\tUnitTest_LIBRARY      = ${UnitTest_LIBRARY}")
    message(STATUS "\tUnitTest_LIBRARIES    = ${UnitTest_LIBRARIES}")
    message(STATUS "")
  endif() 
endif()    


##############################################################################
# OpenMP - http://openmp.org/
#
# comment out set_feature_info per
# https://software.lanl.gov/ascem/trac/ticket/413#comment:1
##############################################################################
option(ENABLE_OpenMP "Build Amanzi executables with OpenMP" OFF)
#set_feature_info(OpenMP
#                 ENABLE_OpenMP
#                 "OpenMP, multi-platform shared-memory parallel programming"
#                 )
if (ENABLE_OpenMP)
  find_package(OpenMP)
  find_package(OpenMP_Fortran)
endif()


##############################################################################
# PETSc - http://www.mcs.anl.gov/petsc
##############################################################################
#option(ENABLE_PETSC "Enable PETSc APIs in the structured mesh" FALSE)
option(ENABLE_PETSC "Enable PETSc APIs in the structured mesh" ON)
if (ENABLE_Structured OR ENABLE_ALQUIMIA OR ENABLE_PETSC) # FIXME: Sloppy.
  find_package(PETSc)
  if (NOT PETSc_FOUND)
    message(WARNING "Failed to locate PETSc")
  else()
    message(STATUS "PETSc Package information")
    message(STATUS "\tPETSc_VERSION      = ${PETSc_VERSION}")
    message(STATUS "\tPETSc_INCLUDE_DIR  = ${PETSc_INCLUDE_DIR}")
    message(STATUS "\tPETSc_INCLUDE_DIRS = ${PETSc_INCLUDE_DIRS}")
    message(STATUS "\tPETSc_LIBRARY      = ${PETSc_LIBRARY}")
    message(STATUS "\tPETSc_LIBRARIES    = ${PETSc_LIBRARIES}")
    print_link_libraries(${PETSc_LIBRARY})
    message(STATUS "")
  endif()
endif()


##############################################################################
# ALQUIMIA / PFLOTRAN 
##############################################################################
if (ENABLE_ALQUIMIA) 
  find_package(PETSc)
  if (NOT PETSC_FOUND) 
    message(WARNING "Failed to locate PETSC")
  endif() 

#  if (ENABLE_PFLOTRAN)
    find_package(PFLOTRAN)
    if (PFLOTRAN_FOUND)
       message(STATUS "PFLOTRAN Package information")
       message(STATUS "\tPFLOTRAN_INCLUDE_DIR  = ${PFLOTRAN_INCLUDE_DIR}")
       message(STATUS "\tPFLOTRAN_INCLUDE_DIRS = ${PFLOTRAN_INCLUDE_DIRS}")
       message(STATUS "\tPFLOTRAN_LIBRARY_DIR  = ${PFLOTRAN_LIBRARY_DIR}")
       message(STATUS "\tPFLOTRAN_LIBRARY      = ${PFLOTRAN_LIBRARY}")
       message(STATUS "\tPFLOTRAN_LIBRARIES    = ${PFLOTRAN_LIBRARIES}")
       print_link_libraries(${PFLOTRAN_LIBRARY})
       message(STATUS "")
    endif()
#  endif()

#  if (ENABLE_CRUNCHTOPE)
     find_package(CrunchTope)
     if (CrunchTope_FOUND)
        message(STATUS "CrunchTope Package information")
        message(STATUS "\tCrunchTope_INCLUDE_DIR  = ${CrunchTope_INCLUDE_DIR}")
        message(STATUS "\tCrunchTope_INCLUDE_DIRS = ${CrunchTope_INCLUDE_DIRS}")
        message(STATUS "\tCrunchTope_LIBRARY_DIR  = ${CrunchTope_LIBRARY_DIR}")
        message(STATUS "\tCrunchTope_LIBRARY      = ${CrunchTope_LIBRARY}")
        message(STATUS "\tCrunchTope_LIBRARIES    = ${CrunchTope_LIBRARIES}")
        print_link_libraries(${CrunchTope_LIBRARY})
        message(STATUS "")
     endif()
#  endif()

  if ((NOT PFLOTRAN_FOUND) AND (NOT CrunchTope_FOUND))
    message(WARNING "Failed to locate either PFLOTRAN or CrunchTope")
  endif()

  find_package(Alquimia)
  if (NOT Alquimia_FOUND)
    message(WARNING "Failed to locate Alquimia")
  else()
    message(STATUS "Alquimia Package information")
    message(STATUS "\tAlquimia_INCLUDE_DIR  = ${Alquimia_INCLUDE_DIR}")
    message(STATUS "\tAlquimia_INCLUDE_DIRS = ${Alquimia_INCLUDE_DIRS}")
    message(STATUS "\tAlquimia_LIBRARY_DIR  = ${Alquimia_LIBRARY_DIR}")
    message(STATUS "\tAlquimia_LIBRARY      = ${Alquimia_LIBRARY}")
    message(STATUS "\tAlquimia_LIBRARIES    = ${Alquimia_LIBRARIES}")
    print_link_libraries(${Alquimia_LIBRARY})
    message(STATUS "")
  endif()
endif()

##############################################################################
# CLM LSM 
##############################################################################
if (ENABLE_CLM)
     find_package(CLM)
     if (CLM_FOUND)
        message(STATUS "CLM Package information")
        message(STATUS "\tCLM_INCLUDE_DIR  = ${CLM_INCLUDE_DIR}")
        message(STATUS "\tCLM_INCLUDE_DIRS = ${CLM_INCLUDE_DIRS}")
        message(STATUS "\tCLM_LIBRARY_DIR  = ${CLM_LIBRARY_DIR}")
        message(STATUS "\tCLM_LIBRARY      = ${CLM_LIBRARY}")
        message(STATUS "\tCLM_LIBRARIES    = ${CLM_LIBRARIES}")
        print_link_libraries(${CLM_LIBRARY})
        message(STATUS "")
     endif()
  endif()
