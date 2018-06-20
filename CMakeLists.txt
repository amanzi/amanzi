# -*- mode: cmake -*-

#
# ATS
#       Root CMakeLists.txt file
#

# Require cmake 3.1.3 or higher for c++11 support
cmake_minimum_required(VERSION 3.1.3)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

# Standard CMake modules
include(CMakeDependentOption)
include(FeatureSummary)


# --------------------------------------------------------------------------- #
# Enable testing
# --------------------------------------------------------------------------- #
enable_testing()
cmake_dependent_option(ENABLE_TESTS "Enable unit testing" ON
                       "ENABLE_UnitTest" ON)
add_feature_info(TESTS
                 ENABLE_TESTS
                 "Toggle for unit tests")
if (ENABLE_TESTS)
  set(BUILD_TESTS 1)
endif()    

# --------------------------------------------------------------------------- #
# Define the project name
# This command will define
#     ATS_SOURCE_DIR
#     ATS_BINARY_DIR
# --------------------------------------------------------------------------- #
project(ATS)

message(STATUS "\n-------------------------------------------")
message(STATUS "\n-- CMake: Configuring ATS build/install.\n--")
message(STATUS "----------------------------------------")


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
  message(STATUS "Amanzi_SOURCE_DIR=${Amanzi_SOURCE_DIR}")
else()
  message(FATAL_ERROR "Failed to locate Amanzi")
endif()

# Useful variables pointing to directories in the source tree
set(AMANZI_SOURCE_TOOLS_DIR    "${AMANZI_SOURCE_DIR}/tools")

# add amanzi cmake
set(AMANZI_MODULE_PATH "${Amanzi_SOURCE_DIR}/tools/cmake")
set(CMAKE_MODULE_PATH 
    ${AMANZI_MODULE_PATH}
    ${AMANZI_MODULE_PATH}/Modules
    ${AMANZI_MODULE_PATH}/Utils)


# --------------------------------------------------------------------------- #
#  Define the compilers
# --------------------------------------------------------------------------- #
set(CMAKE_C_COMPILER ${Amanzi_C_COMPILER})
set(CMAKE_CXX_COMPILER ${Amanzi_CXX_COMPILER})
set(CMAKE_Fortran_COMPILER ${Amanzi_Fortran_COMPILER})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${Amanzi_C_COMPILER_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Amanzi_CXX_COMPILER_FLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${Amanzi_Fortran_COMPILER_FLAGS}")

if (CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions(-DENABLE_DBC)
endif()

message(STATUS " --- CMAKE_C_COMPILER ${CMAKE_C_COMPILER} ")
message(STATUS " --- CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER} ")
message(STATUS " --- CMAKE_Fortran_COMPILER ${CMAKE_Fortran_COMPILER} ")


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


if (ENABLE_ALQUIMIA)
  # Use Alquimia as a chemistry engine.
  add_definitions(-DALQUIMIA_ENABLED)
endif()

message(STATUS "Silo Enabled?: ${Amanzi_TPL_Silo_ENABLED}")
message(STATUS "Silo Enabled?: ${ENABLE_Silo}")
message(STATUS "Alquimia Enabled?: ${Amanzi_TPL_Alquimia_ENABLED}")
message(STATUS "Alquimia Enabled?: ${ENABLE_Alquimia}")

include_directories(${Amanzi_TPL_Trilinos_INCLUDE_DIRS})
link_directories(${Amanzi_TPL_Trilinos_LIBRARY_DIRS})

# add the macros that help with registration of factories
include(${Amanzi_LIBRARY_DIR}/RegisterEvaluators.cmake)


# --------------------------------------------------------------------------- #
#  Testing
# --------------------------------------------------------------------------- #
include(TestManager)


add_subdirectory(src)



