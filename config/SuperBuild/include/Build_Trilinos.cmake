#  -*- mode: cmake -*-

#
# Build TPL: Trilinos
#    
# --- Define all the directories and common external project flags
set(trilinos_depend_projects NetCDF ExodusII Boost)
if(ENABLE_HYPRE)
  list(APPEND trilinos_depend_projects HYPRE)
endif()
define_external_project_args(Trilinos
                             TARGET trilinos
                             DEPENDS ${trilinos_depend_projects})

# --- Define the configuration parameters   

#  - Trilinos Package Configuration

if(Trilinos_Build_Config_File)
  message(STATUS "Including Trilinos build configuration file ${Trilinos_Build_Config_File}")
  if ( NOT EXISTS ${Trilinos_Build_Config_File} )
    message(FATAL_ERROR "File ${Trilinos_Build_Config_File} does not exist.")
  endif()
  include(${Trilinos_Build_Config_File})
endif()

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
    list(APPEND Trilinos_CMAKE_TPL_ARGS "-D${mpi_var}:STRING=${${mpi_var}}")
  endif()
endforeach() 

# BLAS
option(ENABLE_BLA_Search "Search for BLAS/LAPACK libraries" OFF)
if (ENABLE_BLA_Search)
    if ( NOT BLA_VENDOR )
      set(BLA_VENDOR All)
    endif()
    message(STATUS "Searching BLAS libraries vendor - ${BLA_VENDOR}")
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
if (ENABLE_BLA_Search)
    if ( NOT BLA_VENDOR )
      set(BLA_VENDOR All)
    endif()
    message(STATUS "Searching LAPACK libraries vendor - ${BLA_VENDOR}")
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
            message(STATUS "Trilinos LAPACK libraries: ${LAPACK_LIBRARIES}")    
endif()

# Boost
list(APPEND Trilinos_CMAKE_TPL_ARGS
            "-DTPL_ENABLE_Boost:BOOL=ON" 
            "-DBoost_INCLUDE_DIRS:FILEPATH=${TPL_INSTALL_PREFIX}/include"
            "-DBoost_LIBRARY_DIRS:FILEPATH=${TPL_INSTALL_PREFIX}/lib")

# NetCDF
list(APPEND Trilinos_CMAKE_TPL_ARGS
            "-DTPL_ENABLE_Netcdf:BOOL=ON"
            "-DNetcdf_INCLUDE_DIRS:FILEPATH=${TPL_INSTALL_PREFIX}/include"
            "-DNetcdf_LIBRARY_DIRS:FILEPATH=${TPL_INSTALL_PREFIX}/lib")

# ExodusII 
list(APPEND Trilinos_CMAKE_TPL_ARGS
            "-DTPL_ENABLE_ExodusII:BOOL=ON" 
            "-DExodusII_LIBRARY_DIRS:FILEPATH=${TPL_INSTALL_PREFIX}/lib"
            "-DExodusII_INCLUDE_DIRS:FILEPATH=${TPL_INSTALL_PREFIX}/include")

# HYPRE
if( ENABLE_HYPRE )
  list(APPEND Trilinos_CMAKE_TPL_ARGS
              "-DTPL_ENABLE_HYPRE:BOOL=ON"
              "-DHYPRE_LIBRARY_DIRS:FILEPATH=${TPL_INSTALL_PREFIX}/lib"
              "-DHYPRE_INCLUDE_DIRS:FILEPATH=${TPL_INSTALL_PREFIX}/include")
endif()

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
#print_variable(Trilinos_Build_Config_File)
#if(Trilinos_Build_Config_File)
#    set(Trilinos_Config_File_ARGS 
#        "-C${Trilinos_Build_Config_File}")
#    print_variable(Trilinos_Config_File_ARGS)
#    message(STATUS "Will add ${Trilinos_Build_Config_File} to the Trilinos configure")    
#endif()    


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
print_variable(Trilinos_CMAKE_LANG_ARGS)

#  --- Define the Trilinos patch step

# Trilinos needs a patch for GNU versions > 4.6
if ( CMAKE_CXX_COMPILER_VERSION )
  if ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" )
    if ( ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "4.6" )
      set(ENABLE_Trilinos_Patch OFF)
    else()
      message(STATUS "Trilinos requires a patch when using"
                     " GNU ${CMAKE_CXX_COMPILER_VERSION}")
      set(ENABLE_Trilinos_Patch ON)
    endif()
  endif()
endif()  

set(Trilinos_PATCH_COMMAND)
if (ENABLE_Trilinos_Patch)
  configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/trilinos-patch-step.sh.in
               ${Trilinos_prefix_dir}/trilinos-patch-step.sh
               @ONLY)
  set(Trilinos_PATCH_COMMAND sh ${Trilinos_prefix_dir}/trilinos-patch-step.sh)
endif()  
#print_variable(Trilinos_PATCH_COMMAND)

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
		    CMAKE_ARGS          ${Trilnos_Config_File_ARGS}
                    CMAKE_CACHE_ARGS    ${Trilinos_CMAKE_LANG_ARGS} 
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
