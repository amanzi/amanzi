#
# Optional ToolsBuild
#

#
# Set minimum required version of CMake
# NB: Must go before project() call.
#
cmake_minimum_required(VERSION 3.2)

#
# Language standard
# NB: Seems like a good idea to have this before project() call.
#
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ############################################################################ #
# Define project name (ToolsBuild_SOURCE_DIR and ToolsBuild_BINARY_DIR)
# ############################################################################ #
message(STATUS "\n-- CMake: Configuring Amanzi Tools build/install.\n--")
project(ToolsBuild C CXX Fortran)


# ############################################################################ #
# Compiler Definition Warnings
# ############################################################################ #
if (NOT CMAKE_C_COMPILER)
  message(WARNING "C compiler not specified. CMake will guess!")
endif()

if (NOT CMAKE_CXX_COMPILER)
  message(WARNING "C++ compiler not specified. CMake will guess!")
endif()

if (NOT CMAKE_Fortran_COMPILER)
  message(WARNING "Fortran compiler not specified. CMake will guess!")
endif()

#if (NOT SPACK_BINARY)
#  message(WARNING "Spack binary not defined. Disabling spack builds.")
#endif()


# ############################################################################ #
# Define Amanzi directories 
# ############################################################################ #
file(TO_CMAKE_PATH ${ToolsBuild_SOURCE_DIR}/../.. Amanzi_SOURCE_DIR)


# ############################################################################ #
# Update the CMake module path
# ############################################################################ #
set(Amanzi_MODULE_PATH ${Amanzi_SOURCE_DIR}/tools/cmake
                       ${Amanzi_SOURCE_DIR}/tools/cmake/Utils
                       ${Amanzi_SOURCE_DIR}/tools/cmake/Modules)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} 
	              ${SuperBuild_SOURCE_DIR}/cmake
	              ${Amanzi_MODULE_PATH})


# ############################################################################ #
# Location of the Build_* files and template files
# ############################################################################ #
set(ToolsBuild_BUILD_FILES_DIR "${ToolsBuild_SOURCE_DIR}/include")
set(ToolsBuild_TEMPLATE_FILES_DIR "${ToolsBuild_SOURCE_DIR}/templates")


# ############################################################################ #
# Compiler checks
# ############################################################################ #
include(FortranCInterface)
FortranCInterface_VERIFY()

include(DefineCompilerVersion)
define_compiler_version()

include(CheckMPISourceCompiles)
check_mpi_source_compiles(MPI_WRAPPERS_IN_USE)
if (NOT MPI_WRAPPERS_IN_USE)
  # MPI_WRAPPERS_IN_USE is FALSE if any of the compilers fail to compile an MPI app
  #message(FATAL_ERROR "Missing a compiler type for MPI applications -- exiting.
  message(WARING "Missing a compiler type for MPI applications.
                  See previous error messages to determine missing compiler.")
endif()


# ############################################################################ #
# Build options
# ############################################################################ #
set(CMAKE_BUILD_TYPE "Release")
message(STATUS "Setting build type -- ${CMAKE_BUILD_TYPE}")

# Download directory 
if (NOT TOOLS_DOWNLOAD_DIR)
  set(TOOLS_DOWNLOAD_DIR "${ToolsBuild_BINARY_DIR}/Downloads")
endif()
message(STATUS "Download Tools files in ${TOOLS_DOWNLOAD_DIR}") 

# Disable external web downloads
include(${SuperBuild_SOURCE_DIR}/TPLVersions.cmake)
option(DISABLE_EXTERNAL_DOWNLOAD "Flag to disable external web downloads" FALSE)
if (DISABLE_EXTERNAL_DOWNLOAD)
  message(STATUS "Will not download files from the web, instead will search "
                 "for archive files in ${TOOLS_DOWNLOAD_DIR}")
else()
  include(CheckDownload)
  check_download(TEST_FILE MD5SUMS TEST_URL ${AMANZI_TOOLS_DOWNLOAD_URL})
endif()

# Tools install location
if (NOT TOOLS_INSTALL_PREFIX)
  set(TOOLS_INSTALL_PREFIX "${ToolsBuild_BINARY_DIR}/external-projects")
endif()
message(STATUS "Install Tools in ${TOOLS_INSTALL_PREFIX}")

# also set the CMAKE_INSTALL_PREFIX accordingly
set(CMAKE_INSTALL_PREFIX ${TOOLS_INSTALL_PREFIX})


# ############################################################################ #
# MPI Configuration 
# ############################################################################ #
# --- The FindMPI package for CMake for versions 2.8.6> use compiler wrappers
#     to determine include paths, link flags and libraries for each language.
if (MPI_WRAPPERS_IN_USE)
  set(MPI_C_COMPILER ${CMAKE_C_COMPILER})
  set(MPI_CXX_COMPILER ${CMAKE_CXX_COMPILER})
  set(MPI_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})

else()
  message(STATUS "Searching for MPI compiler wrappers to define"
                 " compile and link flags")

  find_program(MPI_C_COMPILER
               NAME mpicc
               HINTS ${MPI_INSTALL_PREFIX} ENV MPIROOT ENV MPI_ROOT ENV MPI_HOME ENV MPIHOME ENV MPICH_ROOT 
               PATH_SUFFIXES bin Bin
               )

  find_program(MPI_CXX_COMPILER
               NAMES mpicxx mpic++ mpiCC
               HINTS ${MPI_INSTALL_PREFIX} ENV MPIROOT ENV MPI_ROOT ENV MPI_HOME ENV MPIHOME ENV MPICH_ROOT 
               PATH_SUFFIXES bin Bin
               )

  find_program(MPI_Fortran_COMPILER
               NAMES mpif90
               HINTS ${MPI_INSTALL_PREFIX} ENV MPIROOT ENV MPI_ROOT ENV MPI_HOME ENV MPIHOME ENV MPICH_ROOT 
               PATH_SUFFIXES bin Bin
               )
endif()

# --- This package find call will set
#     MPI_<lang>_FOUND
#     MPI_<lang>_COMPILE_FLAGS
#     MPI_<lang>_INCLUDE_PATH
#     MPI_<lang>_LINK_FLAGS
#     MPI_<lang>_LIBRARIES
#
#     MPIEXEC
#     MPIEXEC_NUMPROC_FLAG
#     MPIEXEC_PREFLAGS
#     MPIEXEC_POSTFLAGS
find_package(MPI)


# The FindMPI package uses FindPackageHandleStandardArgs
# to MPI_<lang>_FOUND. For consistency, I set MPI_FORTRAN_FOUND
# to MPI_Fortran_FOUND
set(MPI_Fortran_FOUND ${MPI_FORTRAN_FOUND})

# Warn the user if MPI information is not found
if (NOT MPI_C_FOUND)
  message(WARNING "Failed to locate MPI C include and library files")
endif()

if (NOT MPI_CXX_FOUND)
  message(WARNING "Failed to locate MPI C++ include and library files")
endif()

if (NOT MPI_Fortran_FOUND)
  message(WARNING "Failed to locate MPI Fortran include and library files")
endif()

# Add MPI to the build if any MPI_*_FOUND == FALSE
if (MPI_C_FOUND AND MPI_CXX_FOUND AND MPI_Fortran_FOUND)
  message(STATUS "Found MPI")
  set(BUILD_MPI FALSE)
  set(mpi_var_list COMPILER COMPILE_FLAGS INCLUDE_PATH LINK_FLAGS LIBRARIES)
  foreach ( lang C CXX Fortran )
    foreach ( var ${mpi_var_list} )
      set(var MPI_${lang}_${var})
      if ("${var}")
	message(STATUS "      ${var}: ${${var}}")
      endif()
    endforeach()
  endforeach()  

else()
  message(STATUS "Will build MPI with the other Tools")
  set(BUILD_MPI TRUE)
  set(MPI_C_COMPILER        ${TOOLS_INSTALL_PREFIX}/bin/mpicc)
  set(MPI_CXX_COMPILER      ${TOOLS_INSTALL_PREFIX}/bin/mpicxx)
  set(MPI_Fortran_COMPILER  ${TOOLS_INSTALL_PREFIX}/bin/mpif90)
  set(MPIEXEC               ${TOOLS_INSTALL_PREFIX}/bin/mpirun)
  set(MPI_EXEC              ${TOOLS_INSTALL_PREFIX}/bin/mpirun)
endif()

# --- Add the include paths, link flags and libraries if using a 
#     user requested MPI but not the wrappers from this MPI installation.
if (NOT MPI_WRAPPERS_IN_USE AND NOT BUILD_MPI)
  foreach (lang C CXX Fortran)
    # Build MPI compile flags
    set(add_mpi_flags_list)
    if (MPI_${lang}_COMPILE_FLAGS)
      foreach (flag ${MPI_${lang}_COMPILE_FLAGS} )
	list(APPEND add_mpi_flags_list ${flag})
      endforeach()
    endif()

    if (MPI_${lang}_INCLUDE_PATH)
      foreach (dir ${MPI_${lang}_INCLUDE_PATH} )
	list(APPEND add_mpi_flags_list "-I${dir}")
      endforeach()
    endif() 

    if (MPI_${lang}_LINK_FLAGS)
      foreach (flag ${MPI_${lang}_LINK_FLAGS} )
	list(APPEND add_mpi_flags_list ${flag})
      endforeach()
    endif()

    list(REMOVE_DUPLICATES add_mpi_flags_list)
    list(LENGTH add_mpi_flags_list num_flags)

    if ("${num_flags}" GREATER "0")
      include(BuildWhitespaceString)
      if (CMAKE_${lang}_FLAGS)
	build_whitespace_string(new_compile_flags ${CMAKE_${lang}_FLAGS} ${add_mpi_flags_list})
      else()	
	build_whitespace_string(new_compile_flags ${add_mpi_flags_list})
      endif()
      set(CMAKE_${lang}_FLAGS "${new_compile_flags}")
      message(STATUS "Update CMAKE_${lang}_FLAGS: ${CMAKE_${lang}_FLAGS}")
    endif()  
  endforeach()
endif()

# --- Amanzi uses MPI_EXEC* not MPIEXEC* variables. This allows the user to 
#     override the find package results.

# - MPI execute binary
if (NOT MPI_EXEC)
  if (MPIEXEC)
    set(MPI_EXEC "${MPIEXEC}" CACHE STRING "Select MPI executable from FindMPI results")
  else()
    find_program(MPI_EXEC
                 NAMES mpirun mpiexec aprun openmpirun
                 HINTS ENV MPI_ROOT ENV MPIROOT ENV MPI_PREFIX ENV MPI_HOME ENV MPIHOME)
  endif()
endif()

# - Number of MPI ranks flag
set(MPI_EXEC_NUMPROCS_FLAG_DFLT -n)
if (NOT MPI_EXEC_NUMPROCS_FLAG)
  if (MPIEXEC_NUMPROC_FLAG)
    set(MPI_EXEC_NUMPROCS_FLAG "${MPIEXEC_NUMPROC_FLAG}" CACHE STRING "Set MPI number of procs flag from FindMPI")
  else()
    set(MPI_EXEC_NUMPROCS_FLAG ${MPI_EXEC_NUMPROCS_FLAG_DFLT})
  endif()
endif()

# - Maximum number of processors. This is a limit for the test suite
#   Some tests require too many processors and it increases the execution time
#   considerably. 
set(MPI_EXEC_MAX_NUMPROCS_DFLT 8)
if (NOT MPI_EXEC_MAX_NUMPROCS)
  include(ProcessorCount)
  ProcessorCount(proc_count)
  if (NOT proc_count EQUAL 0)
    math(EXPR MPI_EXEC_MAX_NUMPROCS "${proc_count} * 2") 
    message(STATUS "Detected ${proc_count} processors and will set maximum to ${MPI_EXEC_MAX_NUMPROCS}")
  else()
    set(MPI_EXEC_MAX_NUMPROCS ${MPI_EXEC_MAX_NUMPROCS_DFLT})
  endif()
endif()  

# - Set the pre and post flags
#   Usage:
#   ${MPI_EXEC} ${MPI_EXEC_NUMPROCS_FLAG} PROCS ${MPI_EXEC_PREFLAGS} EXECUTABLE ${MPI_EXEC_POSTFLAGS}
if (NOT MPI_EXEC_PREFLAGS)
  if (MPIEXEC_PREFLAGS)
    set(MPI_EXEC_PREFLAGS "${MPIEXEC_PRFLAGS}" CACHE STRING "Set MPI execute pre flags")
  endif()
endif()

if (NOT MPI_EXEC_POSTFLAGS)
  if (MPIEXEC_POSTFLAGS)
    set(MPI_EXEC_POSTFLAGS "${MPIEXEC_PRFLAGS}" CACHE STRING "Set MPI execute post flags")
  endif()
endif()


# ############################################################################ #
# Set common build compiler flags, build types and directories
# ############################################################################ #

# CMake compiler settings for any package built with CMake

# C language flags
set(Amanzi_CMAKE_C_COMPILER_ARGS
     -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
     -DCMAKE_C_FLAGS_DEBUG:STRING=${CMAKE_C_FLAGS_DEBUG}
     -DCMAKE_C_FLAGS_MINSIZEREL:STRING=${CMAKE_C_FLAGS_MINSIZEREL}
     -DCMAKE_C_FLAGS_RELEASE:STRING=${CMAKE_C_FLAGS_RELEASE}
     -DCMAKE_C_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_C_FLAGS_RELWITHDEBINFO})

# C++ language flags 
set(Amanzi_CMAKE_CXX_COMPILER_ARGS
     -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
     -DCMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}
     -DCMAKE_CXX_FLAGS_MINSIZEREL:STRING=${CMAKE_CXX_FLAGS_MINSIZEREL}
     -DCMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}
     -DCMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_CXX_FLAGS_RELWITHDEBINFO})

# Fortran language flags    
set(Amanzi_CMAKE_Fortran_COMPILER_ARGS
     -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS}
     -DCMAKE_Fortran_FLAGS_DEBUG:STRING=${CMAKE_Fortran_FLAGS_DEBUG}
     -DCMAKE_Fortran_FLAGS_MINSIZEREL:STRING=${CMAKE_Fortran_FLAGS_MINSIZEREL}
     -DCMAKE_Fortran_FLAGS_RELEASE:STRING=${CMAKE_Fortran_FLAGS_RELEASE}
     -DCMAKE_Fortran_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_Fortran_FLAGS_RELWITHDEBINFO})

# Link flags
set(Amanzi_CMAKE_LINKER_ARGS
     -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
     -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS})

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
message(STATUS "\tC flags\t\t${Amanzi_COMMON_CFLAGS}")
message(STATUS "\tC++ flags\t${Amanzi_COMMON_CXXFLAGS}")
message(STATUS "\tFortran flags\t${Amanzi_COMMON_FCFLAGS}")

# Common link flags for projects that do not use CMake
set(Amanzi_COMMON_LDFLAGS)
if (CMAKE_EXE_LINK_FLAGS)
  build_whitespace_string(Amanzi_COMMON_LDFLAGS ${CMAKE_EXE_LINK_FLAGS})
endif()


# ############################################################################ #
# Begin MPI builds 
# ############################################################################ #

# --- Add .NOTPARALLEL target in order to force the top level make to 
#     execute serially

add_custom_target(.NOTPARALLEL)

# Create configure, build, install and test targets for each tool
set_property(DIRECTORY PROPERTY
             EP_STEP_TARGETS download patch configure build install test)

# Include the CMake module ExternalProject
include(ExternalProject)

# Include the macro that defines uniform build, install and Add_ExternalProject args
include(DefineExternalProjectArgs)

# MPI
if (BUILD_MPI)
  include(${ToolsBuild_BUILD_FILES_DIR}/Build_OpenMPI.cmake)
  set(MPI_PROJECT OpenMPI)
endif()


# ############################################################################ #
#  Final steps
# ############################################################################ #
#  --- Create a CMake cache script
set(ToolsBuild_Amanzi_CACHE_SCRIPT amanzi-tools-config.cmake)
configure_file(${ToolsBuild_TEMPLATE_FILES_DIR}/${ToolsBuild_Amanzi_CACHE_SCRIPT}.in
               ${ToolsBuild_BINARY_DIR}/${ToolsBuild_Amanzi_CACHE_SCRIPT}
               @ONLY)
install(FILES "${ToolsBuild_BINARY_DIR}/${ToolsBuild_Amanzi_CACHE_SCRIPT}"             
        DESTINATION ${TOOLS_INSTALL_PREFIX}/share/cmake)

