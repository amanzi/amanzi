# -*- mode: cmake -*-
# 
# Amanzi Third Party Library (TPL) Definitions
#

# Standard CMake modules see CMAKE_ROOT/Modules
include(FeatureSummary)

# Amanzi CMake modules see <root source>/tools/cmake
include(CheckMPISourceCompiles)
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
find_package(Trilinos 10.6 REQUIRED
             HINTS ${Trilinos_DIR}
             PATH_SUFFIXES include)
            
if ( Trilinos_FOUND )

    message(STATUS "Found Trilinos: ${Trilinos_LIBRARY_DIR}")

    # Amanzi uses Epetra and Teuchos utils throughout the code. 
    # This find_package call defines Epetra_* variables.
    # Amanzi developers should use these variables
    # for libraries that ONLY use Epetra/Teuchos and avoid
    # using the ALL POWERFUL(TM) Trilinos_LIBRARIES.
    # When/If we create wrappers, using this variable
    # will make that transition easier.
    find_package(Epetra
                 NO_MODULE
                 HINTS ${Trilinos_DIR}
                 PATH_SUFFIXES include
                 )
                
    find_package(Teuchos
                 NO_MODULE
                 HINTS ${Trilinos_DIR}
                 PATH_SUFFIXES include
                 )
                
    # STK (Mesh framework) is not a default Trilinos package
    # Must explicity build this package. This is a check to see
    # if STK is in Trilinos. 
    if (ENABLE_STK_Mesh)
      find_package(STK 
        NO_MODULE
        HINTS ${Trilinos_DIR}
        PATH_SUFFIXES include
        )
    endif()

    # NOX non-linear solver used in flow
    find_package(NOX
                 NO_MODULE
                 HINTS ${Trilinos_DIR}
                 PATH_SUFFIXES include
                 )
                
    # For some reason, Trilinos defines dependent TPLs in *_TPL_LIBRARIES not
    # in *_LIBRARIES. We update the variables so the usage of these variables
    # is consistent with other FindXXX modules.
    list(APPEND Epetra_LIBRARIES "${Epetra_TPL_LIBRARIES}")
    list(APPEND Epetra_INCLUDE_DIRS "${Epetra_TPL_INCLUDE_DIRS}")
    list(APPEND Teuchos_LIBRARIES "${Teuchos_TPL_LIBRARIES}")
    list(APPEND Teuchos_INCLUDE_DIRS "${Teuchos_TPL_INCLUDE_DIRS}")
    list(APPEND STK_LIBRARIES "${STK_TPL_LIBRARIES}")
    list(APPEND STK_INCLUDE_DIRS "${STK_TPL_INCLUDE_DIRS}")
    list(APPEND NOX_LIBRARIES "${NOX_TPL_LIBRARIES}")
    list(APPEND NOX_INCLUDE_DIRS "${NOX_TPL_INCLUDE_DIRS}")

    list(APPEND Trilinos_LIBRARIES "${Trilinos_TPL_LIBRARIES}")
    list(APPEND Trilinos_INCLUDE_DIRS "${Trilinos_TPL_INCLUDE_DIRS}")
else()
    message(FATAL_ERROR "Can not locate Trilinos configuration file\n"
                        " Please define the location of your Trilinos installation\n"
                        "using -D Trilinos_DIR:FILEPATH=<install path>\n")
endif()    

# Trilinos can contain 20 or more libraries (packages). The variable Trilinos_LIBRARIES
# does not have full path names only library names. I suspect if it did include
# full path names the link command would exceed the command line character limit on
# many platforms, thus we need to add Trilinos library path to the link to find these
# libraries. Since Amanzi calls many Trilinos packages directly we'll add Trilinos
# to all link commands here. Yes, this is overkill. We'll have wrappers someday.
link_directories(${Trilinos_LIBRARY_DIRS})

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
