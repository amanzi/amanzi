# -*- mode: cmake -*-
#
# Set up defines necessary to build BoxLib-based code

if(__CCSE_OPTIONS_INCLUDED)
  return()
endif()
set(__CCSE_OPTIONS_INCLUDED 1)

include(FortranCInterface)
include(${FortranCInterface_BINARY_DIR}/Output.cmake)

#message(STATUS "BoxLib-specific compile settings:")
if(FortranCInterface_GLOBAL_SUFFIX STREQUAL ""  AND FortranCInterface_GLOBAL_CASE STREQUAL "UPPER")
#    message(STATUS "   Fortran name mangling scheme to UPPERCASE (upper case, no append underscore)")
    set(BL_FORTLINK UPPERCASE)
elseif(FortranCInterface_GLOBAL_SUFFIX STREQUAL ""  AND FortranCInterface_GLOBAL_CASE STREQUAL "LOWER")
#    message(STATUS "   Fortran name mangling scheme to LOWERCASE (lower case, no append underscore)")
    set(BL_FORTLINK LOWERCASE)
elseif(FortranCInterface_GLOBAL_SUFFIX STREQUAL "_" AND FortranCInterface_GLOBAL_CASE STREQUAL "LOWER")
#    message(STATUS "   Fortran name mangling scheme to UNDERSCORE (lower case, append underscore)")
    set(BL_FORTLINK "UNDERSCORE")
#else()
#    message(AUTHOR_WARNING "Fortran to C mangling not backward compatible with older style BoxLib code") 
endif()

set(BL_MACHINE ${CMAKE_SYSTEM_NAME})

set(BL_DEFINES "BL_NOLINEVALUES;BL_PARALLEL_IO;BL_SPACEDIM=${AMANZI_SPACEDIM};BL_FORT_USE_${BL_FORTLINK};BL_${BL_MACHINE};BL_USE_${AMANZI_PRECISION};MG_USE_FBOXLIB;MG_USE_F90_SOLVERS;MG_USE_FBOXLIB")

if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  list(APPEND BL_DEFINES NDEBUG)
elseif ("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
  list(APPEND BL_DEFINES NDEBUG)
elseif ("${CMAKE_BUILD_TYPE}" STREQUAL "MinSizeRel")
  list(APPEND BL_DEFINES NDEBUG)
endif()

if (ENABLE_MPI)
  # bandre: I think the amanzi config requires that the mpi compilers
  # be set through the CC/CXX/FC environment variables before cmake is
  # called. This is overwriting those values and causing the incorrect
  # values to be used?
  #find_package(MPI REQUIRED)
  list(APPEND BL_DEFINES BL_USE_MPI)
  set(CMAKE_CC_FLAGS "${CMAKE_CC_FLAGS} ${MPI_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_FLAGS}")
  set(CCSE_LIBDIR_MPI_SUFFIX .MPI)
else()
  set(CCSE_LIBDIR_MPI_SUFFIX)
endif()

if (ENABLE_OpenMP)
  set(CCSE_LIBDIR_OMP_SUFFIX .OMP)
else()
  set(CCSE_LIBDIR_OMP_SUFFIX)
endif(ENABLE_OpenMP)

get_directory_property(defs COMPILE_DEFINITIONS)
if(defs)
  list(APPEND BL_DEFINES ${defs})
endif()


if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wno-deprecated")
endif(CMAKE_COMPILER_IS_GNUCXX)

list(APPEND BL_DEFINES "AMANZI")

if (ENABLE_PETSC)

  set(PETSc_DIR $ENV{PETSc_DIR})
  if ("${PETSc_DIR}" STREQUAL "")
    message(FATAL_ERROR "Must define env variable PETSC_DIR if ENABLE_PETSC=ON")
  endif()

  # This adds undesirable path to the top of the search list and leads to conflicts. 
  #if (${APPLE})
  #  include_directories(/usr/local/include) # For Homebrew valgrind installs
  #endif()
  include_directories(${PETSc_DIR}/include)
  list(APPEND BL_DEFINES BL_USE_PETSC)
  set(PETSC_LIB_DIR ${PETSc_DIR}/lib)
  link_directories(${PETSc_LIB_DIR})
  set(PETSC_LIBS petsc)
  set(PETSC_EXT_LIBS X11)
endif()

set_directory_properties(PROPERTIES COMPILE_DEFINITIONS "${BL_DEFINES}")

