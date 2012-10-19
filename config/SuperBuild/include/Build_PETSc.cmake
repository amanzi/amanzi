#  -*- mode: cmake -*-

#
# Build TPL:  PETSc 
#    
# --- Define all the directories and common external project flags
define_external_project_args(PETSc TARGET petsc BUILD_IN_SOURCE)

# --- Download packages PETSc needs
set(petsc_packages ParMetis SuperLU SuperLUDist)
get_filename_component(real_download_path ${TPL_DOWNLOAD_DIR} REALPATH)

message(STATUS "Checking PETSc required packages: ${petsc_packages}")
foreach ( _pack ${petsc_packages} )
  set(_url      "${${_pack}_URL_STRING}")
  set(_archive  "${${_pack}_ARCHIVE_FILE}")
  set(_md5sum   "${${_pack}_MD5_SUM}")
  if ( EXISTS "${real_download_path}/${_archive}" )
    message(STATUS "\tFound ${_archive}")
  else()
    if (DISABLE_EXTERNAL_DOWNLOAD)
      message(FATAL_ERROR "You have disabled external downloads, however"
	                  " ${real_download_path}/${_archive} does not exist")
    else()
      message(STATUS "Downloading ${_archive} for PETSc")
      file(DOWNLOAD  "${_url}/${_archive}" "${real_download_path}/${_archive}"
                     SHOW_PROGRESS
                     INACTIVITY_TIMEOUT 180
                     EXPECTED_MD5SUM ${_md5sum})
    endif()
  endif()  
endforeach()

         

# --- Define configure parameters

# Use the common cflags, cxxflags
include(BuildWhitespaceString)
build_whitespace_string(petsc_cflags
                       ${Amanzi_COMMON_CFLAGS})

build_whitespace_string(petsc_cxxflags
                       ${Amanzi_COMMON_CXXFLAGS})
set(cpp_flag_list 
    ${Amanzi_COMMON_CFLAGS}
    ${Amanzi_COMMON_CXXFLAGS})
list(REMOVE_DUPLICATES cpp_flag_list)
build_whitespace_string(petsc_cppflags ${cpp_flags_list})

# Set PETSc debug flag
if ( "${CMAKE_BUILD_TYPE}" STREQUAL "Release" )
  set(petsc_debug_flag 0)
else()
  set(petsc_debug_flag 1)
endif()

# PETSc MPI flag
set(petsc_mpi_flag --with-mpi=1)

# PETSc SuperLU flags
# For now we allow PETSc to download and build this package
# It should be a separate TPL. Error with the download or
# building of this package will appear to be an error in the
# petsc-configure target. See the log files for more detailed
# information.
set(petsc_superlu_flags 
         --download-superlu_dist=${real_download_path}/${SuperLUDist_ARCHIVE_FILE}
         --download-parmetis=${real_download_path}/${ParMetis_ARCHIVE_FILE}
         --download-superlu=${real_download_path}/${SuperLU_ARCHIVE_FILE})

# PETSc install directory
set(petsc_install_dir ${TPL_INSTALL_PREFIX}/${PETSc_BUILD_TARGET}-${PETSc_VERSION})

if (DEFINED LAPACK_LIBRARIES) 
  set(PETSC_LAPACK_OPTION "--with-lapack-lib=${LAPACK_LIBRARIES}")
else()
  set(PETSC_LAPACK_OPTION )
endif()
if (DEFINED BLAS_LIBRARIES) 
  set(PETSC_BLAS_OPTION "--with-blas-lib=${BLAS_LIBRARIES}")
else()
  set(PETSC_BLAS_OPTION )
endif()


# --- Add external project build 
ExternalProject_Add(${PETSc_BUILD_TARGET}
                    DEPENDS   ${PETSc_PACKAGE_DEPENDS}             # Package dependency target
                    TMP_DIR   ${PETSc_tmp_dir}                     # Temporary files directory
                    STAMP_DIR ${PETSc_stamp_dir}                   # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}               # Download directory
                    URL          ${PETSc_URL}                      # URL may be a web site OR a local file
                    URL_MD5      ${PETSc_MD5_SUM}                  # md5sum of the archive file
                    # -- Configure
                    SOURCE_DIR        ${PETSc_source_dir}          # Source directory
                    CONFIGURE_COMMAND
                              <SOURCE_DIR>/configure
                                          --prefix=<INSTALL_DIR>
                                          --with-cc=${CMAKE_C_COMPILER}
                                          --with-cxx=${CMAKE_CXX_COMPILER}
                                          --with-fc=${CMAKE_Fortran_COMPILER}
                                          --CFLAGS=${petsc_cflags}
                                          --CXXFLAGS=${petsc_cxxflags}
                                          --with-debugging=${petsc_debug_flag}
                                          ${PETSC_LAPACK_OPTION}
                                          ${PETSC_BLAS_OPTION}
                                          ${petsc_mpi_flag}
                                          ${petsc_superlu_flags}
                    # -- Build
                    BINARY_DIR        ${PETSc_build_dir}           # Build directory 
                    BUILD_COMMAND     $(MAKE)                      # Run the CMake script to build
                    BUILD_IN_SOURCE   ${PETSc_BUILD_IN_SOURCE}     # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${petsc_install_dir}  # Install directory, NOT in the usual directory
                    # -- Output control
                    ${PETSc_logging_args})

# --- Useful variables for other packages that depend on PETSc
set(PETSC_DIR ${petsc_install_dir})
