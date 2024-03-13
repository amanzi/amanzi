# -*- mode: cmake -*-
#
#
# Amanzi Build Options
#
#
# This file is intended define build options
# related to compile options, build types, etc.
# Options related to Third Party Libraries (TPL)
# can be found in AmanziTPL.cmake

# Standard CMake modules
include(CMakeDependentOption)
include(FeatureSummary)

enable_language(C)
enable_language(CXX)
enable_language(Fortran)

# No idea why we need this.
# I think it was required for Franklin build. -- lpritch
if (PREFER_STATIC_LIBRARIES)
  # Prefer static libraries, but don't require that everything must be static. 
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a .lib)
endif(PREFER_STATIC_LIBRARIES)

if (BUILD_STATIC_EXECUTABLES)
    set(CMAKE_EXE_LINKER_FLAGS -static)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)       # remove -Wl,-Bdynamic
    set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
    set(CMAKE_SHARED_LIBRARY_C_FLAGS)         # remove -fPIC
    set(CMAKE_SHARED_LIBRARY_CXX_FLAGS)
    set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)    # remove -rdynamic
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
endif(BUILD_STATIC_EXECUTABLES)

#
# Options
# 

# DBC - Design by contract
option(ENABLE_DBC "Enable Design By Contract (DBC) checking" ON)
add_feature_info(DBC
                 ENABLE_DBC
                 "Toggle design by contract (DBC) checking")
if (ENABLE_DBC)
  add_definitions("-DENABLE_DBC")
endif()    

# Trilinos error checking is defined in macros
# if ( ${CMAKE_BUILD_TYPE} STREQUAL "Debug" )
#     add_definitions("-DHAVE_FATAL_MESSAGES:BOOL=TRUE")
# endif()


# Testing
# We do not have a consistent way to activate the unit and other tests
# should have a single switch for this. -- lpritch
#cmake_dependent_option(ENABLE_TESTS "Enable unit testing" ON
#                       "ENABLE_UnitTest" ON)
if(ENABLE_TESTS OR ENABLE_UnitTest) 
  set(ENABLE_UnitTest ON)
  set(ENABLE_TEST ON)
endif()

add_feature_info(TESTS
                 ENABLE_TESTS
                 "Toggle for unit tests")
if (ENABLE_TESTS)
  set(BUILD_TESTS 1)
endif()    

if (DEFINED ENABLE_PETSC)
else()
  set(ENABLE_PETSC FALSE CACHE BOOL "Enables the PETSc library." FORCE)
endif()

# Some platforms require all binaries linking to MPI
# only run through the MPIEXEC binary
option(TESTS_REQUIRE_MPIEXEC "Run all tests with the MPIEXEC binary" FALSE)

# Need this option if the PATH environment does not include '.'
option(TESTS_REQUIRE_FULLPATH "Append full path to test binaries" TRUE)

# This is needed to allow the consistent registration of evaluators with
# their factory, see src/common/standalone_simulation_coordinator/CMakeLists.txt,
# src/common/standalone_simulation_coordinator/Main.cc, 
# src/state/CMakeLists.txt, and src/0state/cell_volume_evaluator_reg.hh
# for how to use this
include(RegisterEvaluators)

# Two sets of maintained time integrators
option(ENABLE_NewTimeIntegrators "Build Amanzi with dev time integrators" FALSE)

# Turn off physics, just build the framework
option(ENABLE_AmanziPhysicsModule "Build Amanzi with native physics libraries" TRUE)
option(ENABLE_ATSPhysicsModule "Build Amanzi with ATS physics libraries" FALSE)


