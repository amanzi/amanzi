#  -*- mode: cmake -*-

#
# Build TPL:  PETSc 
#    
# --- Define all the directories and common external project flags
define_external_project_args(PETSc TARGET petsc BUILD_IN_SOURCE)

# --- Define configure parameters

# Use the common cflags, cxxflags
include(BuildWhitespaceString)
build_whitespace_string(petsc_cflags
                       -I${TPL_INSTALL_PREFIX}/include
                       ${Amanzi_COMMON_CFLAGS})

build_whitespace_string(petsc_cxxflags
                       -I${TPL_INSTALL_PREFIX}/include
                       ${Amanzi_COMMON_CXXFLAGS})
set(cpp_flag_list 
    -I${TPL_INSTALL_PREFIX}/include
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
         --download-superlu_dist
	 --download-parmetis
	 --download-superlu)

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
                                          --prefix=${TPL_INSTALL_PREFIX}
                                          --with-cc=${CMAKE_C_COMPILER}
                                          --with-cxx=${CMAKE_CXX_COMPILER}
                                          --with-fc=${CMAKE_Fortran_COMPILER}
                                          --CFLAGS=${petsc_cflags}
                                          --CXXFLAGS=${petsc_cxxflags}
                                          --with-debugging=${petsc_debug_flag}
                                          ${petsc_mpi_flag}
					  ${petsc_superlu_flags}
                    # -- Build
                    BINARY_DIR        ${PETSc_build_dir}           # Build directory 
                    BUILD_COMMAND     $(MAKE)                      # Run the CMake script to build
                    BUILD_IN_SOURCE   ${PETSc_BUILD_IN_SOURCE}     # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}         # Install directory
                    # -- Output control
                    ${PETSc_logging_args})
