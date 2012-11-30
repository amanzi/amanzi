# -*- mode: cmake -*-
# 
# Amanzi Third Party Library (TPL) Definitions
#

# Standard CMake modules see CMAKE_ROOT/Modules
include(FeatureSummary)

# Amanzi CMake modules see <root source>/tools/cmake
include(CheckMPISourceCompiles)
include(TrilinosMacros)
include(PrintVariable)


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
find_package( Boost COMPONENTS system filesystem program_options regex REQUIRED)
set_feature_info(Boost
                 "C++ Extension library"
                 "http://www.boost.org"
                 "Required by the MPC")

if ( Boost_VERSION) 

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

# We need to use the project-local HDF5 finder. Temporarily Change
# policy CMP0017 if were using cmake 2.8.3 or later
if (${ADJUST_POLICY})
  cmake_policy(SET CMP0017 OLD)
endif()

find_package(HDF5 1.8.0 REQUIRED)
if ( NOT HDF5_IS_PARALLEL ) 
    message(WARNING     "The HDF5 installation found in ${HDF5_DIR} is not "
                        "a parallel build. At this time, this installation "
                        "is compatible with other TPLs. Soon Amanzi will "
                        "require a parallel enabled HDF5. Please update your "
                        "HDF5 installation to include MPI I/O symbols"
            )            
endif(NOT HDF5_IS_PARALLEL)
set_feature_info(HDF5
                "I/O library that creates HDF5 formatted files"
                "http://www.hdfgroup.org/HDF5"
                "Required library for several components in Amanzi"
                )

# Restore policy of preferring offical CMake modules over local ones.
if (${ADJUST_POLICY})
  cmake_policy(SET CMP0017 NEW)
endif()

##############################################################################
# Trilinos http://trilinos.sandia.gov
##############################################################################
# This command alters Trilinos_DIR. If it finds the configuration file
# Trilinos_DIR is set to the path the configuration file was found.
if ( NOT Trilinos_INSTALL_PREFIX )
  message(WARNING "Use Trilinos_INSTALL_PREFIX"
                  " to define the Trilinos installation location"
		  "\n-DTrilinos_INSTALL_PREFIX:PATH=<trilnos directory>\n")
endif()
set(Trilinos_MINIMUM_VERSION 11.0.3)
find_package(Trilinos ${Trilinos_MINIMUM_VERSION} REQUIRED
             PATHS ${Trilinos_INSTALL_PREFIX}
             PATH_SUFFIXES include)
            
if ( Trilinos_FOUND )

    message(STATUS "Found Trilinos: ${Trilinos_DIR} (${Trilinos_VERSION})")
    trilinos_package_enabled_tpls(Trilinos)           

    if ( "${Trilinos_VERSION}" VERSION_LESS ${Trilinos_MINIMUM_VERSION} ) 
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
    if ( ENABLE_Unstructured )
      list(APPEND Trilinos_REQUIRED_PACKAGE_LIST NOX ML)
    endif()

    foreach(tri_package ${Trilinos_REQUIRED_PACKAGE_LIST})
      find_package(${tri_package} REQUIRED
                   NO_MODULE
                   HINTS ${Trilinos_INSTALL_PREFIX}
                   PATH_SUFFIXES include lib)
      trilinos_package_enabled_tpls(${tri_package})
      message(STATUS "Located Trilinos package ${tri_package}: ${${tri_package}_DIR}")
      # Update the <PACKAGE>_INCLUDE_DIRS variable 
      foreach( _inc ${${tri_package}_TPL_INCLUDE_DIRS})
	list(APPEND ${tri_package}_INCLUDE_DIRS "${_inc}")
      endforeach()

    endforeach()

    # Now check optional Trilinos packages

    # STK - mesh 
    if ( ENABLE_STK_Mesh )
      find_package(STK
                   NO_MODULE
                   HINTS ${Trilinos_INSTALL_PREFIX}
                   PATH_SUFFIXES include lib)
      if (STK_FOUND)
	message(STATUS "Located Trilinos package STK: ${STK_DIR}")
        trilinos_package_enabled_tpls(STK)
	foreach( _inc "${STK_TPL_INCLUDE_DIRS}")
	    list(APPEND STK_INCLUDE_DIRS "${_inc}")
	endforeach()
      else()  
        message(WARNING "Could not locate STK in ${Trilinos_DIR}. Will not enable STK_Mesh")
        set(ENABLE_STK_Mesh OFF CACHE BOOL "Disable STK Mesh capability" FORCE)
      endif()

    endif()

    # Zoltan - required by MSTK mesh class 
    if ( ENABLE_MSTK_Mesh )
      find_package(Zoltan
                   NO_MODULE
                   HINTS ${Trilinos_INSTALL_PREFIX}
                   PATH_SUFFIXES include lib)
      if (Zoltan_FOUND)
	message(STATUS "Located Trilinos package Zoltan: ${Zoltan_DIR}")
        trilinos_package_enabled_tpls(Zoltan)
	foreach( _inc "${ZOLTAN_TPL_INCLUDE_DIRS}")
	    list(APPEND ZOLTAN_INCLUDE_DIRS "${_inc}")
	endforeach()
      else()  
        message(WARNING "Could not locate Zoltan in ${Trilinos_DIR}. Will not enable MSTK_Mesh")
        set(ENABLE_MSTK_Mesh OFF CACHE BOOL "Disable MSTK Mesh capability" FORCE)
      endif()

    endif()


    option(ENABLE_HYPRE "Enable HYPRE APIs in flow" ON)
    if ( ENABLE_HYPRE )

      # Ifpack - preconditioner package that serves as a wrapper for HYPRE
      find_package(Ifpack 
                   NO_MODULE
                   HINTS ${Trilinos_INSTALL_PREFIX}
                   PATH_SUFFIXES include lib
                   )
      if (Ifpack_FOUND)
	message(STATUS "Located Triilnos package Ifpack: ${Ifpack_DIR}")
        trilinos_package_enabled_tpls(Ifpack)
	foreach( _inc "${Ifpack_TPL_INCLUDE_DIRS}")
	    list(APPEND Ifpack_INCLUDE_DIRS "${_inc}")
	endforeach()
      else()
        message(SEND_ERROR "Trilinos in ${Trilinos_DIR} does not have the Ifpack package")
      endif()

      if ( NOT Ifpack_ENABLE_HYPRE )
        message(WARNING "ENABLE_HYPRE requires the Trilinos package Ifpack with enabled HYPRE."
                           " Deactivating HYPRE APIs")
        set(ENABLE_HYPRE OFF CACHE BOOL "Disable the HYPRE APIs" FORCE)                 
      endif()
    endif()
       
    # Now update the Trilinos_LIBRARIES and INCLUDE_DIRS
    foreach( _inc "${Trilinos_TPL_INCLUDE_DIRS}")
      list(APPEND Trilinos_INCLUDE_DIRS "${_inc}")
    endforeach()

else()
    message(FATAL_ERROR "Can not locate Trilinos configuration file\n"
                        " Please define the location of your Trilinos installation\n"
                        "using -D Trilinos_DIR:FILEPATH=<install path>\n")
endif()    

##############################################################################
# NetCDF - http://www.unidata.ucar.edu/software/netcdf/
##############################################################################
find_package(NetCDF REQUIRED)
set_feature_info(NetCDF
                 "Network Common Data Format (NetCDF)"
                 "http://www.unidata.ucar.edu/software/netcdf/"
                 "Required by ExodusII library")


##############################################################################
# Exodus II -http://sourceforge.net/projects/exodusii
##############################################################################
find_package(ExodusII REQUIRED)
set_feature_info(ExodusII
                 "File format library. Originated from Sandia."
                 "http://sourceforge.net/projects/exodusii/"
                 "Required by all the mesh frameworks to read mesh files")


##############################################################################
# CCSE - http://ccse.lbl.gov/Software/ccse_core.html
##############################################################################
if (ENABLE_Structured)
  find_package(CCSE REQUIRED)
  set_feature_info(CCSE
                   "CCSE BoxLib softare library required for structured grid")
endif()

##############################################################################
############################ Option Processing ###############################
##############################################################################







##############################################################################
#---------------------------- Mesh Frameworks -------------------------------#
##############################################################################

# Enable ALL possible mesh frameworks
#option(ENABLE_ALL_Mesh "Build all Amanzi mesh frameworks" OFF)
#if(ENABLE_ALL_Mesh)
#    set(ENABLE_STK_Mesh ON)
#    set(ENABLE_MOAB_Mesh ON)
#    set(ENABLE_MSTK_Mesh ON)
#endif()    
#set_feature_info(ALL_Mesh
#                 ENABLE_ALL_Mesh
#                 "Build all available mesh frameworks"
#                  )    

##############################################################################
# STK - Sierra Mesh Tool Kit part of Trilinos
##############################################################################
option(ENABLE_STK_Mesh  "Build Amanzi with the STK mesh framework" OFF)
set_feature_info(STK_Mesh
                 ENABLE_STK_Mesh
                 "Sierra Mesh Tool Kit (STK Mesh) a Trilinos package"
                 )


##############################################################################
# MOAB - svn co https://svn.mcs.anl.gov/repos/ITAPS/MOAB/trunk MOAB
##############################################################################
option(ENABLE_MOAB_Mesh "Build Amanzi with the MOAB mesh framework" OFF)
set_feature_info(MOAB_Mesh
                 ENABLE_MOAB_Mesh
                 "A Mesh-Oriented datABase"
                 )
if (ENABLE_MOAB_Mesh)
    find_package(MOAB REQUIRED)
endif()

##############################################################################
# MSTK - https://software.lanl.gov/MeshTools/trac/raw-attachment/wiki/WikiStart/mstk-1.80.tar.gz
##############################################################################
option(ENABLE_MSTK_Mesh "Build Amanzi with the MOAB mesh framework" OFF)
set_feature_info(MSTK_Mesh
                 ENABLE_MSTK_Mesh
                 "A mesh framework"
                 )
if (ENABLE_MSTK_Mesh)
    find_package(MSTK REQUIRED)
endif() 





##############################################################################
#-------------------------- Optional Libraries ------------------------------#
##############################################################################

##############################################################################
# ASCEMIO - http://www.cgns.sourceforge.net/
##############################################################################
option(ENABLE_ASCEMIO  "Build Amanzi output library with ASCEM-IO parallelIO" OFF)
set_feature_info(ASCEMIO
                  ENABLE_ASCEMIO
                 "ASCEM-IO Scalable Parallel I/O module for Environmental Management Applications"
                 "http://ascem-io.secure-water.org"
                 "Required to produce VisIt files in parallel"
                 )
#if (ENABLE_ASCEMIO)
if (ENABLE_Unstructured)
    find_package(ASCEMIO REQUIRED)
else()
    find_package(ASCEMIO)
endif() 

##############################################################################
# CGNS - http://www.cgns.sourceforge.net/
##############################################################################
option(ENABLE_CGNS  "Build Amanzi output library with CGNS" OFF)
set_feature_info(CGNS
                  ENABLE_CGNS
                 "CFD General Notation System"
                 "http://cgns.sourceforge.net"
                 "Required to produce VisIt files"
                 )
if (ENABLE_CGNS)
    find_package(CGNS REQUIRED)
endif() 

##############################################################################
# UnitTest++ - http://unittest-cpp.sourceforge.net/
##############################################################################
option(ENABLE_UnitTest "Build Amanzi unit tests. Requires UnitTest++" ON)
set_feature_info(UnitTest
                 ENABLE_UnitTest
                 "C++ unit test framework"
                 )
if (ENABLE_UnitTest)
    find_package(UnitTest)
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
endif()

##############################################################################
# PETSc - http://www.mcs.anl.gov/petsc
##############################################################################

option(ENABLE_PETSC "Enable PETSc APIs in the structured mesh" FALSE)
if ( ENABLE_Structured AND  ENABLE_PETSC )
  find_package(PETSc)
  if ( NOT PETSC_FOUND )
    message(WARNING "Failed to locate PETSc")
  else()
    message(STATUS "PETSc Package information")
    message(STATUS "\tPETSC_VERSION       =${PETSC_VERSION}")
    message(STATUS "\tPETSC_INCLUDE_DIR   =${PETSC_INCLUDE_DIR}")
    message(STATUS "\tPETSC_INCLUDE_DIRS  =${PETSC_INCLUDE_DIRS}")
    message(STATUS "\tPETSC_LIBRARY       =${PETSC_LIBRARY}")
    message(STATUS "\tPETSC_LIBRARIES     =${PETSC_LIBRARIES}")
  endif()
endif()


