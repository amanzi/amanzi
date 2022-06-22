#
# Optional ToolsBuild
#

#
# Set minimum required version of CMake
# NB: Must go before project() call.
#
cmake_minimum_required(VERSION 3.10.0)

#
# Language standard
# NB: Seems like a good idea to have this before project() call.
#
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ############################################################################ #
# Define project name (ToolsBuild_SOURCE_DIR and ToolsBuild_BINARY_DIR)
# ############################################################################ #
message(STATUS "------------------------------------------------")
message(STATUS "\n-- CMake: Configuring Amanzi Tools build/install\n--")
message(STATUS "------------------------------------------------")
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
	              ${ToolsBuild_SOURCE_DIR}/cmake
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
include(${ToolsBuild_SOURCE_DIR}/../SuperBuild/TPLVersions.cmake)
option(DISABLE_EXTERNAL_DOWNLOAD "Flag to disable external web downloads" FALSE)
if (DISABLE_EXTERNAL_DOWNLOAD)
  message(STATUS "Will not download files from the web, instead will search "
                 "for archive files in ${TOOLS_DOWNLOAD_DIR}")
else()
  include(CheckDownload)
  check_download(TEST_FILE MD5SUMS TEST_URL ${AMANZI_TPLS_DOWNLOAD_URL})
endif()

# Tools install location
if (NOT TOOLS_INSTALL_PREFIX)
  set(TOOLS_INSTALL_PREFIX "${ToolsBuild_BINARY_DIR}/external-projects")
endif()
message(STATUS "Install Tools in ${TOOLS_INSTALL_PREFIX}")

# also set the CMAKE_INSTALL_PREFIX accordingly
set(CMAKE_INSTALL_PREFIX ${TOOLS_INSTALL_PREFIX})


# ############################################################################ #
# Set common build compiler flags, build types and directories
# ############################################################################ #

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

set(MPI_C_COMPILER        ${TOOLS_INSTALL_PREFIX}/bin/mpicc)
set(MPI_CXX_COMPILER      ${TOOLS_INSTALL_PREFIX}/bin/mpicxx)
set(MPI_Fortran_COMPILER  ${TOOLS_INSTALL_PREFIX}/bin/mpif90)
set(MPIEXEC               ${TOOLS_INSTALL_PREFIX}/bin/mpirun)
set(MPI_EXEC              ${TOOLS_INSTALL_PREFIX}/bin/mpirun)

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
if (${TOOLS_MPI} STREQUAL "mpich") 
  include(${ToolsBuild_BUILD_FILES_DIR}/Build_MPICH.cmake)
else()
  include(${ToolsBuild_BUILD_FILES_DIR}/Build_OpenMPI.cmake)
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

