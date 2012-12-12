#  -*- mode: cmake -*-

#
# Build TPL:  HYPRE 
#    
# --- Define all the directories and common external project flags
define_external_project_args(HYPRE 
                             TARGET hypre
                             DEPENDS ${MPI_PROJECT}
                             BUILD_IN_SOURCE)

# --- Define configure parameters

# Use the common cflags, cxxflags
include(BuildWhitespaceString)
build_whitespace_string(hypre_cflags
                       -I${TPL_INSTALL_PREFIX}/include
                       ${Amanzi_COMMON_CFLAGS})

build_whitespace_string(hypre_cxxflags
                       -I${TPL_INSTALL_PREFIX}/include
                       ${Amanzi_COMMON_CXXFLAGS})
set(cpp_flag_list 
    -I${TPL_INSTALL_PREFIX}/include
    ${Amanzi_COMMON_CFLAGS}
    ${Amanzi_COMMON_CXXFLAGS})
list(REMOVE_DUPLICATES cpp_flag_list)
build_whitespace_string(hypre_cppflags ${cpp_flags_list})

# Disable OpenMP with HYPRE for now
# Is OpenMP available
#if ( ENABLE_OpenMP )
# find_package(OpenMP)

 set(hypre_openmp_opt)
# if ( OPENMP_FOUND )
#   set(hypre_openmp_opt --with-openmp)
# endif()
#else()
 set(hypre_openmp_opt --without-openmp)
#endif()
print_variable(hypre_openmp_opt)

# Locate LAPACK and BLAS

set(hypre_blas_opt)
find_package(BLAS)
if ( BLAS_FOUND )
  set(hypre_blas_opt --with-blas)
endif()

set(hypre_lapack_opt)
find_package(LAPACK)
if ( LAPACK_FOUND )
  set(hypre_lapack_opt --with-lapack)
endif()

set(hypre_fortran_opt --disable-fortran)

# Build the configure script
set(HYPRE_sh_configure ${HYPRE_prefix_dir}/hypre-configure-step.sh)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/hypre-configure-step.sh.in
               ${HYPRE_sh_configure}
	       @ONLY)

# Configure the CMake command file
set(HYPRE_cmake_configure ${HYPRE_prefix_dir}/hypre-configure-step.cmake)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/hypre-configure-step.cmake.in
               ${HYPRE_cmake_configure}
	       @ONLY)
set(HYPRE_CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${HYPRE_cmake_configure})	

# --- Define the build command

# Build the build script
set(HYPRE_sh_build ${HYPRE_prefix_dir}/hypre-build-step.sh)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/hypre-build-step.sh.in
               ${HYPRE_sh_build}
	       @ONLY)

# Configure the CMake command file
set(HYPRE_cmake_build ${HYPRE_prefix_dir}/hypre-build-step.cmake)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/hypre-build-step.cmake.in
               ${HYPRE_cmake_build}
	       @ONLY)
set(HYPRE_BUILD_COMMAND ${CMAKE_COMMAND} -P ${HYPRE_cmake_build})	

# --- Define the install command

# Build the install script
set(HYPRE_sh_install ${HYPRE_prefix_dir}/hypre-install-step.sh)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/hypre-install-step.sh.in
               ${HYPRE_sh_install}
	       @ONLY)

# Configure the CMake command file
set(HYPRE_cmake_install ${HYPRE_prefix_dir}/hypre-install-step.cmake)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/hypre-install-step.cmake.in
               ${HYPRE_cmake_install}
	       @ONLY)
set(HYPRE_INSTALL_COMMAND ${CMAKE_COMMAND} -P ${HYPRE_cmake_install})	

# --- Add external project build and tie to the ZLIB build target
ExternalProject_Add(${HYPRE_BUILD_TARGET}
                    DEPENDS   ${HYPRE_PACKAGE_DEPENDS}             # Package dependency target
                    TMP_DIR   ${HYPRE_tmp_dir}                     # Temporary files directory
                    STAMP_DIR ${HYPRE_stamp_dir}                   # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}               # Download directory
                    URL          ${HYPRE_URL}                      # URL may be a web site OR a local file
                    URL_MD5      ${HYPRE_MD5_SUM}                  # md5sum of the archive file
                    # -- Configure
                    SOURCE_DIR        ${HYPRE_source_dir}          # Source directory
                    CONFIGURE_COMMAND ${HYPRE_CONFIGURE_COMMAND}
                    # -- Build
                    BINARY_DIR        ${HYPRE_build_dir}           # Build directory 
                    BUILD_COMMAND     ${HYPRE_BUILD_COMMAND}       # Run the CMake script to build
                    BUILD_IN_SOURCE   ${HYPRE_BUILD_IN_SOURCE}     # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${TPL_INSTALL_PREFIX}         # Install directory
		    INSTALL_COMMAND  ${HYPRE_INSTALL_COMMAND}
                    # -- Output control
                    ${HYPRE_logging_args})

# --- Useful variables that depend on HYPRE
include(BuildLibraryName)
set(HYPRE_INCLUDE_DIRS "${TPL_INSTALL_PREFIX}/include")
build_library_name(HYPRE HYPRE_LIBRARY APPEND_PATH ${TPL_INSTALL_PREFIX}/lib)
set(HYPRE_LIBRARIES    "${HYPRE_LIBRARY}")
