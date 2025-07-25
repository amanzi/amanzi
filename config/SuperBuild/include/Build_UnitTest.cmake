#  -*- mode: cmake -*-

#
# Build TPL: UnitTest
# 
# --- Define all the directories and common external project flags
if ( NOT ENABLE_XSDK )
    set(unittest_depend_projects ZLIB)
else()
    set(unittest_depend_projects XSDK)
endif()

define_external_project_args(UnitTest
                             TARGET unittest
                             DEPENDS ${unittest_depend_projects})

# add version version to the autogenerated tpl_versions.h file
amanzi_tpl_version_write(FILENAME ${TPL_VERSIONS_INCLUDE_FILE}
  PREFIX UnitTest
  VERSION ${UnitTest_VERSION_MAJOR} ${UnitTest_VERSION_MINOR} ${UnitTest_VERSION_PATCH})

# --- define the configuration parameters

set(Unittest_CMAKE_PACKAGE_ARGS "")

set(Unittest_CMAKE_TPL_ARGS)

# Pass the following MPI arguments to unittest
set(MPI_CMAKE_ARGS DIR EXEC EXEC_NUMPROCS_FLAG EXE_MAX_NUMPROCS C_COMPILER)
foreach (var ${MPI_CMAKE_ARGS} )
  set(mpi_var "MPI_${var}")
  if ( ${mpi_var} )
    list(APPEND Unittest_CMAKE_TPL_ARGS "-D${mpi_var}:STRING=${${mpi_var}}")
  endif()
endforeach() 

# build type
if ( CMAKE_BUILD_TYPE )
  list(APPEND Unittest_CMAKE_EXTRA_ARGS
              "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}")
endif()

#  - Add CMake configuration file
if(Unittest_Build_Config_File)
    list(APPEND Unittest_Config_File_ARGS
        "-C${Unittest_Build_Config_File}")
    message(STATUS "Will add ${Unittest_Build_Config_File} to the Unittest configure")    
    message(DEBUG "Unittest_CMAKE_EXTRA_ARGS = ${Unittest_CMAKE_EXTRA_ARGS}")
endif()    

# --- Patch the original code
set(UnitTest_patch_file unittest-cmake.patch unittest-testrunner.patch)
set(UnitTest_sh_patch ${UnitTest_prefix_dir}/unittest-patch-step.sh)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/unittest-patch-step.sh.in
               ${UnitTest_sh_patch}
               @ONLY)
# configure the CMake patch step
set(UnitTest_cmake_patch ${UnitTest_prefix_dir}/unittest-patch-step.cmake)
configure_file(${SuperBuild_TEMPLATE_FILES_DIR}/unittest-patch-step.cmake.in
               ${UnitTest_cmake_patch}
               @ONLY)
# set the patch command
set(UnitTest_PATCH_COMMAND ${CMAKE_COMMAND} -P ${UnitTest_cmake_patch})

# --- Add external project build 
set(Unittest_CMAKE_ARGS 
   ${Unittest_CMAKE_PACKAGE_ARGS}
   ${Unittest_CMAKE_TPL_ARGS}
   ${Unittest_CMAKE_EXTRA_ARGS})

# --- Override minimum version
# --- Cache Args may need the type to be float? 
if(CMAKE_MAJOR_VERSION VERSION_EQUAL "4")
  list(APPEND Unittest_CMAKE_ARGS "-DCMAKE_POLICY_VERSION_MINIMUM:STRING=3.5")
endif()
 
ExternalProject_add(${UnitTest_BUILD_TARGET}
                    DEPENDS   ${UnitTest_PACKAGE_DEPENDS}          # Package dependency target
                    TMP_DIR   ${UnitTest_tmp_dir}                  # Temporary files directory
                    STAMP_DIR ${UnitTest_stamp_dir}                # Timestamp and log directory
                    # -- Download and URL definitions
                    DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}
                    URL          ${UnitTest_URL}                   # URL may be a web site OR a local file
                    URL_MD5      ${UnitTest_MD5_SUM}               # md5sum of the archive file
                    PATCH_COMMAND ${UnitTest_PATCH_COMMAND}        # Mods to source
                    # -- Configure
                    SOURCE_DIR       ${UnitTest_source_dir}        # Defining forces CMake to mkdir SOURCE_DIR
                    CMAKE_ARGS       ${Unittest_Config_File_ARGS}
                    CMAKE_CACHE_ARGS ${AMANZI_CMAKE_CACHE_ARGS}    # Ensure uniform build
                                     ${Unittest_CMAKE_ARGS}
                                     -DCMAKE_INSTALL_PREFIX:PATH=${TPL_INSTALL_PREFIX}
                                     -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
                                     -DCMAKE_C_FLAGS:STRING=${Amanzi_COMMON_CFLAGS}  # Ensure uniform build
                                     -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                                     -DCMAKE_CXX_FLAGS:STRING=${Amanzi_COMMON_CXXFLAGS}
                                     -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                                     -DCMAKE_Fortran_FLAGS:STRING=${Amanzi_COMMON_FCFLAGS}
                                     -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
                    # -- Build
                    BUILD_COMMAND     $(MAKE)                       # Run make in build directory $(MAKE) enables parallel build
                    BINARY_DIR        ${UnitTest_build_dir}         # Define the build directory
                    BUILD_IN_SOURCE   ${UnitTest_BUILD_IN_SOURCE}   # Flag in/out source build
                    # -- Install
                    INSTALL_DIR       ${TPL_INSTALL_PREFIX}         # Install directory
                    # -- Output control
                    ${UnitTest_logging_args})

