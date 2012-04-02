#  -*- mode: cmake -*-

#
# Build TPL: Trilinos
#    
# --- Define all the directories and common external project flags
define_external_project_args(Trilinos
                             TARGET trilinos
                             DEPENDS NetCDF ExodusII Boost)

# --- Define the configuration parameters   

#  - Trilinos Package Configuration

# List of packages enabled in the Trilinos build
set(Trilinos_PACKAGE_LIST Teuchos Epetra NOX)
if ( ENABLE_STK_Mesh )
  list(APPEND Trilinos_PACKAGE_LIST STK)
endif()

# Add the Seacas package for Trilinos versions greater than 10.8
if ( "${Trilinos_VERSION}" VERSION_GREATER 10.8 )
  list(APPEND Trilinos_PACKAGE_LIST SEACAS)
endif()

# Generate the Trilinos Package CMake Arguments
set(Trilinos_CMAKE_PACKAGE_ARGS "-DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF")
foreach(package ${Trilinos_PACKAGE_LIST})
  list(APPEND Trilinos_CMAKE_PACKAGE_ARGS "-DTrilinos_ENABLE_${package}:STRING=ON")
endforeach()

# We have had trouble with Trilinos builds that do not set this variable, although
# it seems to build nearly every available package in Trilinos.
list(APPEND Trilinos_CMAKE_PACKAGE_ARGS
            "-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON")

#  - Trilinos TPL Configuration

set(Trilinos_CMAKE_TPL_ARGS)

# MPI
list(APPEND Trilinos_CMAKE_TPL_ARGS "-DTPL_ENABLE_MPI:BOOL=ON")

# Pass the following MPI arguments to Trilinos if they are set 
set(MPI_CMAKE_ARGS DIR EXEC EXEC_NUMPROCS_FLAG EXE_MAX_NUMPROCS C_COMPILER)
foreach (var ${MPI_CMAKE_ARGS} )
  set(mpi_var "MPI_${var}")
  if ( ${mpi_var} )
    list(APPEND Trilinos_CMAKE_TPL_ARGS "-D${mpi_var}=${${mpi_var}}")
  endif()
endforeach() 

# BLAS
option(ENABLE_BLAS_Search "Search for BLAS libraries" OFF)
if (ENABLE_BLAS_Search)
    if ( NOT BLA_VENDOR )
      message(STATUS "Search all possible BLAS vendor types")
    endif()
    find_package(BLAS)
    if (NOT BLAS_FOUND )
      message(FATAL_ERROR "Failed to locate BLAS libraries."
                          "Define BLAS libraries with"
  "\n-DBLAS_LIBRARIES:STRING=....\n"
  "and re-run cmake")
    endif()      
endif()

if ( BLAS_LIBRARIES )
  list(APPEND Trilinos_CMAKE_TPL_ARGS
              "-DTPL_ENABLE_BLAS:BOOL=TRUE")
  list(APPEND Trilinos_CMAKE_TPL_ARGS
              "-DTPL_BLAS_LIBRARIES:FILEPATH=${BLAS_LIBRARIES}")
  message(STATUS "Trilinos BLAS libraries: ${BLAS_LIBRARIES}")    
endif()            
 
# LAPACK
if (ENABLE_BLAS_Search)
    if ( NOT BLA_VENDOR )
      message(STATUS "Search all possible BLAS vendor types")
    endif()
    find_package(LAPACK)
    if ( NOT LAPACK_FOUND )
      message(FATAL_ERROR "Failed to locate LAPACK libraries."
                          "Define LAPACK libraries with"
  "\n-DLAPACK_LIBRARIES:STRING=....\n"
  "and re-run cmake")
    endif()
endif()

if ( LAPACK_LIBRARIES )
  list(APPEND Trilinos_CMAKE_TPL_ARGS
    "-DTPL_LAPACK_LIBRARIES:FILEPATH=${LAPACK_LIBRARIES}")
endif()

# Boost
list(APPEND Trilinos_CMAKE_TPL_ARGS
            "-DTPL_ENABLE_Boost:BOOL=ON" 
            "-DBoost_INCLUDE_DIRS:FILEPATH=${TPL_INSTALL_PREFIX}/include"
            "-DBoost_LIBRARY_DIRS:FILEPATH=${TPL_INSTALL_PREFIX}/lib")

# NetCDF
list(APPEND Trilinos_CMAKE_TPL_ARGS
            "-DTPL_ENABLE_Netcdf:BOOL=ON"
            "-DNetcdf_INCLUDE_DIRS=${TPL_INSTALL_PREFIX}/include"
            "-DNetcdf_LIBRARY_DIRS:PATH=${TPL_INSTALL_PREFIX}/lib")

# ExodusII 
list(APPEND Trilinos_CMAKE_TPL_ARGS
            "-DTPL_ENABLE_ExodusII:BOOL=ON" 
            "-DExodusII_LIBRARY_DIRS:PATH=${TPL_INSTALL_PREFIX}/lib"
            "-DExodusII_INCLUDE_DIRS:PATH=${TPL_INSTALL_PREFIX}/include")

#  - Addtional Trilinos CMake Arguments
set(Trilinos_CMAKE_EXTRA_ARGS
    "-DTrilinos_VERBOSE_CONFIGURE:BOOL=ON"
    "-DTrilinos_ENABLE_TESTS:BOOL=OFF"
    "-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON"
    "-DTrilinos_ENABLE_Fortran:BOOL=OFF"
    )
if ( CMAKE_BUILD_TYPE )
  list(APPEND Trilinos_CMAKE_EXTRA_ARGS
              "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}")
endif()

if ( BUILD_SHARED_LIBS )
  list(APPEND Trilinos_CMAKE_EXTRA_ARGS
    "-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}")
endif()


#  - Add CMake configuration file
if(Trilinos_Build_Config_File)
    list(APPEND Trilinos_CMAKE_EXTRA_ARGS
        "-C${Trilinos_Build_Config_File}")
    message(STATUS "Will add ${Trilinos_Build_Config_File} to the Trilinos configure")    
endif()    


#  - Final Trilinos CMake Arguments 
set(Trilinos_CMAKE_ARGS 
   ${Trilinos_CMAKE_PACKAGE_ARGS}
   ${Trilinos_CMAKE_TPL_ARGS}
   ${Trilinos_CMAKE_EXTRA_ARGS}
   )

# - Final language ARGS
set(Trilinos_CMAKE_LANG_ARGS
                   ${Amanzi_CMAKE_C_COMPILER_ARGS}
                   ${Amanzi_CMAKE_CXX_COMPILER_ARGS})

#  --- Define the Trilinos patch step

# Trilinos needs a patch for GNU versions > 4.6
set(Trilinos_PATCH_COMMAND)
option(ENABLE_Trilinos_Patch "Enable the patch step for the Trilinos build" OFF)
if ( CMAKE_CXX_COMPILER_VERSION )
  if (     ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" )
       AND (${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "4.6" ) )
     message(STATUS "Trilinos requires a patch when using"
                    " GNU ${CMAKE_CXX_COMPILER_VERSION}")
     set(ENABLE_Trilinos_Patch ON)
  endif()
endif()  

if (ENABLE_Trilinos_Patch)
  configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/trilinos-patch-step.sh.in
               ${Trilinos_prefix_dir}/trilinos-patch-step.sh
               @ONLY)
  set(Trilinos_PATCH_COMMAND sh ${Trilinos_prefix_dir}/trilinos-patch-step.sh)
endif()  

# --- Add external project build and tie to the Trilinos build target
ExternalProject_Add(${Trilinos_BUILD_TARGET}
                    DEPENDS   ${Trilinos_PACKAGE_DEPENDS}             # Package dependency target
                    TMP_DIR   ${Trilinos_tmp_dir}                     # Temporary files directory
                    STAMP_DIR ${Trilinos_stamp_dir}                   # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}                  # Download directory
                    URL          ${Trilinos_URL}                      # URL may be a web site OR a local file
                    URL_MD5      ${Trilinos_MD5_SUM}                  # md5sum of the archive file
                    # -- Patch
                    PATCH_COMMAND ${Trilinos_PATCH_COMMAND}
                    # -- Configure
                    SOURCE_DIR    ${Trilinos_source_dir}           # Source directory
                    CMAKE_ARGS    ${Trilinos_CMAKE_LANG_ARGS} 
                                  ${Trilinos_CMAKE_ARGS}
                                 -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
                    # -- Build
                    BINARY_DIR        ${Trilinos_build_dir}        # Build directory 
                    BUILD_COMMAND     $(MAKE)                      # $(MAKE) enables parallel builds through make
                    BUILD_IN_SOURCE   ${Trilinos_BUILD_IN_SOURCE}  # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}        # Install directory
                    # -- Output control
                    ${Trilinos_logging_args})
