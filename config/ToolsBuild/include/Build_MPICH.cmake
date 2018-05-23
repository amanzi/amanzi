#  -*- mode: cmake -*-

#
# Build TPL: MPICH 
# 
# --- Define all the directories and common external project flags
define_external_project_args(MPICH TARGET mpich)

# Build compiler *FLAGS strings. Pick up the CMAKE_BUILD_TYPE flags
include(BuildWhitespaceString)
build_whitespace_string(mpich_cflags ${Amanzi_COMMON_CFLAGS})
build_whitespace_string(mpich_cxxflags ${Amanzi_COMMON_CXXFLAGS})
build_whitespace_string(mpich_fcflags ${Amanzi_COMMON_FCFLAGS})

# --- Add RPATH to the link flags for the compiler wrappers
set(mpich_extra_ldflags "-Wl,-rpath,${TOOLS_INSTALL_PREFIX}/lib")
message(STATUS "mpich_extra_ldflags = ${mpich_extra_ldflags}")
find_package(Threads)

if (APPLE)
  set(mpich_extra_options "--enable-two-level-namespace")
  message(STATUS "MPICH extra options for Darwin: ${mpich_extra_options}")
endif()

# --- Add external project build and tie to the OpenMPI build target
ExternalProject_Add(${MPICH_BUILD_TARGET}
                    DEPENDS   ${MPICH_PACKAGE_DEPENDS}     # Package dependency target
                    TMP_DIR   ${MPICH_tmp_dir}             # Temporary files directory
                    STAMP_DIR ${MPICH_stamp_dir}           # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TOOLS_DOWNLOAD_DIR}
                    URL          ${MPICH_URL}              # URL may be a web site OR a local file
                    URL_MD5      ${MPICH_MD5_SUM}          # md5sum of the archive file
                    # -- Configure
                    SOURCE_DIR   ${MPICH_source_dir}
                    CONFIGURE_COMMAND
                                   <SOURCE_DIR>/configure
                                                --prefix=<INSTALL_DIR>
                                                --enable-fortran
						--enable-shared
						--enable-static
	  					--with-wrapper-ldflags=${mpich_extra_ldflags}
                                                ${mpich_extra_options}
                                                CC=${CMAKE_C_COMPILER}
                                                CXX=${CMAKE_CXX_COMPILER}
                                                FC=${CMAKE_Fortran_COMPILER}
                    # -- Build
                    BINARY_DIR        ${MPICH_build_dir}        # Build directory 
                    BUILD_COMMAND     make -j ${TOOLS_PARALLEL_JOBS}  # $(MAKE) enables parallel builds through make
                    BUILD_IN_SOURCE   ${MPICH_BUILD_IN_SOURCE}  # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TOOLS_INSTALL_PREFIX}
                    # -- Output control
                    ${MPICH_logging_args})

# --- Define variables pointing to compiler wrappers and parallel run commond
set(MPI_C_COMPILER        ${TPL_INSTALL_PREFIX}/bin/mpicc)
set(MPI_CXX_COMPILER      ${TPL_INSTALL_PREFIX}/bin/mpicxx)
set(MPI_Fortran_COMPILER  ${TPL_INSTALL_PREFIX}/bin/mpif90)
set(MPIEXEC               ${TPL_INSTALL_PREFIX}/bin/mpirun)
set(MPI_EXEC              ${TPL_INSTALL_PREFIX}/bin/mpirun)

