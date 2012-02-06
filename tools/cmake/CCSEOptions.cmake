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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS} -g")
if (ENABLE_MPI)
  # bandre: I think the amanzi config requires that the mpi compilers
  # be set through the CC/CXX/FC environment variables before cmake is
  # called. This is overwriting those values and causing the incorrect
  # values to be used?
  #find_package(MPI REQUIRED)
  list(APPEND BL_DEFINES BL_USE_MPI)
  set(CMAKE_CC_FLAGS "${CMAKE_CC_FLAGS} ${MPI_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")
  set(CCSE_LIBDIR_MPI_SUFFIX .MPI)
else()
  set(CCSE_LIBDIR_MPI_SUFFIX)
endif()

if (ENABLE_OpenMP)
  list(APPEND BL_DEFINES BL_USE_OMP)
  find_package(OpenMP)
  set(CMAKE_CC_FLAGS "${CMAKE_CC_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CCSE_LIBDIR_OMP_SUFFIX .OMP)
else()
  set(CCSE_LIBDIR_OMP_SUFFIX)
endif(ENABLE_OpenMP)

get_directory_property(defs COMPILE_DEFINITIONS)
if(defs)
  set_directory_properties(PROPERTIES COMPILE_DEFINITIONS $defs ${BL_DEFINES})
endif()


if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth-64 -Wno-deprecated")
endif(CMAKE_COMPILER_IS_GNUCXX)

if (AMANZI_CHEMEVOL_PKG STREQUAL "AMANZI")
  list(APPEND BL_DEFINES "AMANZI")
elseif (AMANZI_CHEMEVOL_PKG STREQUAL "COREREACT")
  list(APPEND BL_DEFINES "COREREACT")
endif()

set_directory_properties(PROPERTIES COMPILE_DEFINITIONS "${BL_DEFINES}")


