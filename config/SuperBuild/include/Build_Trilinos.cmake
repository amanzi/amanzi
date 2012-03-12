#  -*- mode: cmake -*-

#
# Build TPL: Trilinos
#    
define_external_project(Trilinos 
                        TARGET trilinos
                        DEPENDS netcdf exodusii boost
                        )

# ############################################################################ #
# Trilinos Package Configuration                                               #
# ############################################################################ #

# List of packages enabled in the Trilinos build
set(Trilinos_PACKAGE_LIST Teuchos Epetra NOX STK)

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

# ############################################################################ #
# Trilinos TPL Configuration                                                   #
# ############################################################################ #

set(Trilinos_CMAKE_TPL_ARGS)

# MPI
list(APPEND Trilinos_CMAKE_TPL_ARGS "-DTPL_ENABLE_MPI:BOOL=ON")

# Pass the following MPI arguments to Trilinos if they are set 
set(MPI_CMAKE_ARGS DIR EXEC EXEC_NUMPROCS_FLAG EXE_MAX_NUMPROCS)
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
            "-DBoost_INCLUDE_DIRS:FILEPATH=${Boost_install_dir}/include"
            "-DBoost_LIBRARY_DIRS:FILEPATH=${Boost_install_dir}/lib")

# netCDF
list(APPEND Trilinos_CMAKE_TPL_ARGS
            "-DTPL_ENABLE_Netcdf:BOOL=ON"
            "-DNetcdf_INCLUDE_DIRS=${TPL_INSTALL_PREFIX}/include"
            "-DNetcdf_LIBRARY_DIRS:PATH=${TPL_INSTALL_PREFIX}/lib")

# ExodusII 
list(APPEND Trilinos_CMAKE_TPL_ARGS
            "-DTPL_ENABLE_ExodusII:BOOL=ON" 
            "-DExodusII_LIBRARY_DIRS:PATH=${TPL_INSTALL_PREFIX}/lib"
            "-DExodusII_INCLUDE_DIRS:PATH=${TPL_INSTALL_PREFIX}/include")

# ############################################################################ #
# Addtional Trilinos CMake Arguments                                           #
# ############################################################################ #
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


# ############################################################################ #
# Addtional Trilinos CMake Arguments                                           #
# ############################################################################ #
if(Trilinos_Build_Config_File)
    list(APPEND Trilinos_CMAKE_EXTRA_ARGS
        "-C${Trilinos_Build_Config_File}")
    message(STATUS "Will add ${Trilinos_Build_Config_File} to the Trilinos configure")    
endif()    


# ############################################################################ #
# Trilinos CMake Arguments                                                     #
# ############################################################################ #
set(Trilinos_CMAKE_ARGS 
   ${Trilinos_CMAKE_PACKAGE_ARGS}
   ${Trilinos_CMAKE_TPL_ARGS}
   ${Trilinos_CMAKE_EXTRA_ARGS}
   )

# ############################################################################ #
# Add Trilinos patch                                                           #
# ############################################################################ #

# Trilinos needs a patch for GNU versions > 4.6
# Will leave the option available until this is tested
set(Trilinos_patch_cmd)
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
  configure_file(${SuperBuild_BUILD_FILES_DIR}/trilinos-patch-step.sh.in
               ${Trilinos_prefix_dir}/trilinos-patch-step.sh
               @ONLY)
  set(Trilinos_patch_cmd sh ${Trilinos_prefix_dir}/trilinos-patch-step.sh)
endif()  

# ############################################################################ #
# Add Trilinos as an External Package                                          #
# ############################################################################ #
ExternalProject_Add(${Trilinos_target}
                    DEPENDS boost netcdf exodusii
                    ${Trilinos_ep_directory_args}
                    ${Trilinos_url_args}
                    PATCH_COMMAND ${Trilinos_patch_cmd}
                    CMAKE_ARGS
                              ${Amanzi_CMAKE_COMPILER_ARGS}
                              ${Trilinos_CMAKE_ARGS}
                              -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
                    BUILD_COMMAND $(MAKE)          
                    ${Trilinos_logging_opts}
)
