#  -*- mode: cmake -*-

#
# Build TPL: Trilinos
#    
# --- Define all the directories and common external project flags
set(trilinos_depend_projects ${MPI_PROJECT} NetCDF ExodusII Boost)
if(ENABLE_HYPRE)
  list(APPEND trilinos_depend_projects HYPRE)
endif()
define_external_project_args(Trilinos
                             TARGET trilinos
                             DEPENDS ${trilinos_depend_projects})

# --- Define the configuration parameters   

#  - Trilinos Package Configuration

#if(Trilinos_Build_Config_File)
#  message(STATUS "Including Trilinos build configuration file ${Trilinos_Build_Config_File}")
#  if ( NOT EXISTS ${Trilinos_Build_Config_File} )
#    message(FATAL_ERROR "File ${Trilinos_Build_Config_File} does not exist.")
#  endif()
#  include(${Trilinos_Build_Config_File})
#endif()

# List of packages enabled in the Trilinos build
set(Trilinos_PACKAGE_LIST Teuchos Epetra NOX)
if ( ENABLE_STK_Mesh )
  list(APPEND Trilinos_PACKAGE_LIST STK)
endif()
if ( ENABLE_MSTK_Mesh )
  list(APPEND Trilinos_PACKAGE_LIST Zoltan)
endif()


# Generate the Trilinos Package CMake Arguments
set(Trilinos_CMAKE_PACKAGE_ARGS "-DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF")
foreach(package ${Trilinos_PACKAGE_LIST})
  list(APPEND Trilinos_CMAKE_PACKAGE_ARGS "-DTrilinos_ENABLE_${package}:STRING=ON")
endforeach()

# Remove SEACAS from the build and force STK to use external Exodus
if ( ENABLE_STK_Mesh )
  list(APPEND Trilinos_CMAKE_PACKAGE_ARGS "-DTrilinos_ENABLE_SEACAS:STRING=OFF")
  list(APPEND Trilinos_CMAKE_PACKAGE_ARGS "-DSTK_ENABLE_SEACASExodus:STRING=OFF")
  list(APPEND Trilinos_CMAKE_PACKAGE_ARGS "-DSTK_ENABLE_SEACASNemesis:STRING=OFF")
endif()

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
if ( BLAS_LIBRARIES )
  list(APPEND Trilinos_CMAKE_TPL_ARGS
              "-DTPL_ENABLE_BLAS:BOOL=TRUE")
  list(APPEND Trilinos_CMAKE_TPL_ARGS
              "-DTPL_BLAS_LIBRARIES:STRING=${BLAS_LIBRARIES}")
  message(STATUS "Trilinos BLAS libraries: ${BLAS_LIBRARIES}")    
else()
  message(WARNING "BLAS libraies not set. Trilinos will perform search.") 
endif()            
 
# LAPACK
if ( LAPACK_LIBRARIES )
  list(APPEND Trilinos_CMAKE_TPL_ARGS
              "-DTPL_LAPACK_LIBRARIES:STRING=${LAPACK_LIBRARIES}")
            message(STATUS "Trilinos LAPACK libraries: ${LAPACK_LIBRARIES}")    
else()
  message(WARNING "LAPACK libraies not set. Trilinos will perform search.") 
endif()

# Boost
list(APPEND Trilinos_CMAKE_TPL_ARGS
            "-DTPL_ENABLE_Boost:BOOL=ON" 
            "-DBoost_INCLUDE_DIRS:FILEPATH=${TPL_INSTALL_PREFIX}/include"
            "-DBoost_LIBRARY_DIRS:FILEPATH=${TPL_INSTALL_PREFIX}/lib")

# NetCDF
list(APPEND Trilinos_CMAKE_TPL_ARGS
            "-DTPL_ENABLE_Netcdf:BOOL=ON"
            "-DTPL_Netcdf_INCLUDE_DIRS:STRING=${NetCDF_INCLUDE_DIRS}"
            "-DTPL_Netcdf_LIBRARIES:STRING=${NetCDF_C_LIBRARIES}")


# HYPRE
if( ENABLE_HYPRE )
  list(APPEND Trilinos_CMAKE_TPL_ARGS
              "-DTPL_ENABLE_HYPRE:BOOL=ON"
              "-DTPL_HYPRE_LIBRARIES:STRING=${HYPRE_LIBRARIES}"
              "-DHYPRE_INCLUDE_DIRS:PATH=${TPL_INSTALL_PREFIX}/include")
endif()

#  - Addtional Trilinos CMake Arguments
set(Trilinos_CMAKE_EXTRA_ARGS
    "-DTrilinos_VERBOSE_CONFIGURE:BOOL=ON"
    "-DTrilinos_ENABLE_TESTS:BOOL=OFF"
    "-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON"
    )
if ( CMAKE_BUILD_TYPE )
  list(APPEND Trilinos_CMAKE_EXTRA_ARGS
              "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}")
  message(DEBUG "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
  message(DEBUG "Trilinos_CMAKE_EXTRA_ARGS = ${Trilinos_CMAKE_EXTRA_ARGS}")
endif()

if ( BUILD_SHARED_LIBS )
  list(APPEND Trilinos_CMAKE_EXTRA_ARGS
    "-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}")
  message(DEBUG "Trilinos_CMAKE_EXTRA_ARGS = ${Trilinos_CMAKE_EXTRA_ARGS}")
endif()


#  - Add CMake configuration file
if(Trilinos_Build_Config_File)
    list(APPEND Trilinos_Config_File_ARGS
        "-C${Trilinos_Build_Config_File}")
    message(STATUS "Will add ${Trilinos_Build_Config_File} to the Trilinos configure")    
    message(DEBUG "Trilinos_CMAKE_EXTRA_ARGS = ${Trilinos_CMAKE_EXTRA_ARGS}")
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
		           -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER_USE}
                   ${Amanzi_CMAKE_CXX_COMPILER_ARGS}
		           -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER_USE}
                   ${Amanzi_CMAKE_Fortran_COMPILER_ARGS}
                   -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER_USE})

#  --- Define the Trilinos patch step

# Trilinos needs a patch for GNU versions > 4.6
#LPRITCHif ( CMAKE_CXX_COMPILER_VERSION )
#LPRITCH  if ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" )
#LPRITCH    if ( ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "4.6" )
#LPRITCH      set(ENABLE_Trilinos_Patch OFF)
#LPRITCH    else()
#LPRITCH      message(STATUS "Trilinos requires a patch when using"
#LPRITCH                     " GNU ${CMAKE_CXX_COMPILER_VERSION}")
#LPRITCH      set(ENABLE_Trilinos_Patch ON)
#LPRITCH    endif()
#LPRITCH  endif()
#LPRITCHendif()  
#LPRITCH
#LPRITCHset(Trilinos_PATCH_COMMAND)
#LPRITCHif (ENABLE_Trilinos_Patch)
#LPRITCH    set(Trilinos_patch_file)
#LPRITCH    # Set the patch file name
#LPRITCH    if(CMAKE_CXX_COMPILER_VERSION)
#LPRITCH      if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#LPRITCH        if ( "${CMAKE_CXX_COMPILER_VERSION}" VERSION_LESS "4.6" )
#LPRITCH          message(FATAL_ERROR "ENABLE_Trilinos_Patch is ON, however no patch file exists"
#LPRITCH                              " for version ${CMAKE_CXX_COMPILER_VERSION}.")
#LPRITCH        elseif( "${CMAKE_CXX_COMPILER_VERSION}" VERSION_LESS "4.7" )
#LPRITCH          set(Trilinos_patch_file trilinos-${Trilinos_VERSION}-gcc46.patch)
#LPRITCH        elseif ( "${CMAKE_CXX_COMPILER_VERSION}" VERSION_LESS "4.8" )
#LPRITCH          set(Trilinos_patch_file trilinos-${Trilinos_VERSION}-gcc47.patch)
#LPRITCH        else()
#LPRITCH          message(FATAL_ERROR "ENABLE_Trilinos_Patch is ON, however no patch file exists"
#LPRITCH                             " for version ${CMAKE_CXX_COMPILER_VERSION}.")
#LPRITCH        endif()
#LPRITCH      endif()
#LPRITCH    endif()
#LPRITCH
#LPRITCH    #print_variable(Trilinos_patch_file)
#LPRITCH    if(Trilinos_patch_file)
#LPRITCH       configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/trilinos-patch-step.sh.in
#LPRITCH                      ${Trilinos_prefix_dir}/trilinos-patch-step.sh
#LPRITCH                      @ONLY)
#LPRITCH       set(Trilinos_PATCH_COMMAND sh ${Trilinos_prefix_dir}/trilinos-patch-step.sh)
#LPRITCH    else()
#LPRITCH       message(WARNING "ENABLE_Trilinos_Patch is ON but no patch file found for "
#LPRITCH	               "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} "
#LPRITCH		       "Will not patch Trilinos.")
#LPRITCH    endif()		   
#LPRITCH   		   
#LPRITCHendif()  
#print_variable(Trilinos_PATCH_COMMAND)

# --- Define the Trilinos location
set(Trilinos_install_dir ${TPL_INSTALL_PREFIX}/${Trilinos_BUILD_TARGET}-${Trilinos_VERSION})

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
		    CMAKE_ARGS          ${Trilinos_Config_File_ARGS}
                    CMAKE_CACHE_ARGS    ${Trilinos_CMAKE_LANG_ARGS} 
                                        ${Trilinos_CMAKE_ARGS}
                                        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
					-DTrilinos_ENABLE_Stratimikos:BOOL=FALSE
                    # -- Build
                    BINARY_DIR        ${Trilinos_build_dir}        # Build directory 
                    BUILD_COMMAND     $(MAKE)                      # $(MAKE) enables parallel builds through make
                    BUILD_IN_SOURCE   ${Trilinos_BUILD_IN_SOURCE}  # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${Trilinos_install_dir}        # Install directory
                    # -- Output control
                    ${Trilinos_logging_args})

# --- Useful variables for packages that depends on Trilinos
set(Trilinos_INSTALL_PREFIX  ${Trilinos_install_dir})
set(Zoltan_INSTALL_PREFIX "${Trilinos_install_dir}")
