# -*- mode: cmake -*-

#
# ATS
#       Root CMakeLists.txt file
#

# Require cmake 2.8 or higher
cmake_minimum_required(VERSION 2.8.2)

# Enable testing, ctest needs this
# all add_test commands are ignored unless this is called!
set(BUILD_TESTS true)

# Define the project name
# This command will define
#     ATS_SOURCE_DIR
#     ATS_BINARY_DIR
project(ATS)
message(STATUS "in ats_dir: ATS_SOURCE_DIR=${ATS_SOURCE_DIR}")

# --------------------------------------------------------------------------- #
# Search for Amanzi
#  Set Amanzi_DIR when invoking cmake
# --------------------------------------------------------------------------- #
find_package(Amanzi)

if ( Amanzi_FOUND )
  message(STATUS "Located Amanzi")
  message(STATUS "Amanzi_VERSION=${Amanzi_VERSION}")
  message(STATUS "Amanzi_DIR=${Amanzi_DIR}")
  message(STATUS "Amanzi_INCLUDE_DIR=${Amanzi_INCLUDE_DIR}")      
  message(STATUS "Amanzi_ENABLE_STRUCTURED=${Amanzi_ENABLE_STRUCTURED}")
  message(STATUS "Amanzi_ENABLE_UNSTRUCTURED=${Amanzi_ENABLE_UNSTRUCTURED}")
else()
  message(FATAL_ERROR "Failed to locate Amanzi")
endif()

# --------------------------------------------------------------------------- #
#  Define the compilers
# --------------------------------------------------------------------------- #
#set(CMAKE_C_COMPILER ${Amanzi_C_COMPILER})
#set(CMAKE_CXX_COMPILER ${Amanzi_CXX_COMPILER})
if (${PROFILE_FLAGS})
  set(CMAKE_BUILD_TYPE Release)
  add_definitions("-g -O3")
endif()

if (CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions(-DENABLE_DBC)
endif()


# --------------------------------------------------------------------------- #
#  Define the include paths
# --------------------------------------------------------------------------- #

# Amanzi install include path
include_directories(${Amanzi_INCLUDE_DIR})

# Amanzi TPLs
# The pattern for each variable is
# Amanzi_TPL_<tpl name>_<variable>
# Example: For the Teuchos package in Trilinos
# Amanzi_TPL_Teuchos_DIR          => Location of the TeuchosConfig.cmake file
# Amanzi_TPL_Teuchos_INCLUDE_DIR  => Include path for the Teuchos header files
# Amanzi_TPL_Teuchos_INCLUDE_DIRS => Include paths for Teuchos AND TPLs that Teuchos depends on
# Amanzi_TPL_Teuchos_LIBRARY_DIRS => Library paths for Teuchos library and TPL libraries Teuchos calls
# See AmanziConfigTPL.cmake for the entire list

include_directories(${Amanzi_TPL_Teuchos_INCLUDE_DIRS})
include_directories(${Amanzi_TPL_Epetra_INCLUDE_DIRS})
include_directories(${Amanzi_TPL_ASCEMIO_INCLUDE_DIRS})

# Each TPL is flagged as enabled ON(TRUE)|OFF(FALSE)

if ( Amanzi_TPL_CGNS_ENABLED )
  include_directories(${Amanzi_TPL_CGNS_INCLUDE_DIRS})
  add_definitions(-DENABLE_CGNS)
endif()

if ( Amanzi_TPL_STK_ENABLED )
  include_directories(${Amanzi_TPL_STK_INCLUDE_DIRS})
  add_definitions(-DENABLE_STK_Mesh)
endif()


# enable structured/unstructured
if ( Amanzi_ENABLE_UNSTRUCTURED )
  add_definitions(-DENABLE_Unstructured)
endif()

if ( Amanzi_ENABLE_STRUCTURED )
  add_definitions(-DENABLE_Structured)
endif()

# We make every effort in the Amanzi build system to include TPL libraries 
# in the target definitions. Avoid using Amanzi_TPL_*_LIBRARY or
# Amanzi_TPL_*_LIBRARIES variables. Using these variables when it
# is not neccessary may result in a command line length that exceeds
# the OS limits and will prevent the host code from compiling.
# Unfortunately Trilinos does not provide exported targets
# The host code must explicitly link against the Trilinos library
# path to link against the Trilinos libraries.
link_directories(${Amanzi_TPL_Trilinos_LIBRARY_DIRS})


# --------------------------------------------------------------------------- #
#  Building against Amanzi
# --------------------------------------------------------------------------- #

# Source files for all binaries and libraries found under src
add_subdirectory(src)

