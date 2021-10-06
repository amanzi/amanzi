#  -*- mode: cmake -*-

#
# Build TPL: Ginkgo
#    
define_external_project_args(Ginkgo
                            TARGET ginkgo
                            DEPENDS ${ginkgo_depend_projects})

set(Ginkgo_CMAKE_ARGS "-DGINKGO_BUILD_HWLOC:BOOL=OFF")
list(APPEND Ginkgo_CMAKE_ARGS "-DGINKGO_BUILD_HIP:BOOL=OFF")
list(APPEND Ginkgo_CMAKE_ARGS "-DGINKGO_BUILD_OMP:BOOL=OFF")

if(ENABLE_CUDA)
    list(APPEND Ginkgo_CMAKE_ARGS "-DGINKGO_BUILD_CUDA:BOOL=ON")
endif()

# --- Define the Ginkgo location
set(Ginkgo_install_dir ${TPL_INSTALL_PREFIX})

# --- If downloads are disabled point to local repository
if ( DISABLE_EXTERNAL_DOWNLOAD )
  STRING(REGEX REPLACE ".*\/" "" Ginkgo_GIT_REPOSITORY_LOCAL_DIR ${Ginkgo_GIT_REPOSITORY})
  set (Ginkgo_GIT_REPOSITORY_TEMP ${TPL_DOWNLOAD_DIR}/${Ginkgo_GIT_REPOSITORY_LOCAL_DIR})
else()
  set (HYPRE_GIT_REPOSITORY_TEMP ${Ginkgo_GIT_REPOSITORY})
endif()
message(STATUS "Ginkgo git repository = ${Ginkgo_GIT_REPOSITORY_TEMP}")

message(STATUS "Ginkgo Args = ${Ginkgo_CMAKE_ARGS}")

# --- Add external project build and tie to the Ginkgo build target
ExternalProject_Add(${Ginkgo_BUILD_TARGET}
                    DEPENDS   ${Ginkgo_PACKAGE_DEPENDS}             # Package dependency target
                    TMP_DIR   ${Ginkgo_tmp_dir}                     # Temporary files directory
                    STAMP_DIR ${Ginkgo_stamp_dir}                   # Timestamp and log directory
                    # -- Download and URL definitions
                    GIT_REPOSITORY ${Ginkgo_GIT_REPOSITORY_TEMP}              
                    GIT_TAG        ${Ginkgo_GIT_TAG}      
                    # -- Update (one way to skip this step is use null command)
                    UPDATE_COMMAND ""
                    # -- Patch
                    PATCH_COMMAND ${Ginkgo_PATCH_COMMAND}
                    # -- Configure
                    SOURCE_DIR        ${Ginkgo_source_dir}           # Source directory
                    CMAKE_ARGS        ${Ginkgo_Config_File_ARGS}
                    CMAKE_CACHE_ARGS  ${AMANZI_CMAKE_CACHE_ARGS}   # Ensure uniform build
                                      ${Ginkgo_CMAKE_ARGS} 
                                      #-DGINKGO_BUILD_HIP:BOOL=OFF
                                      #-DGINKGO_BUILD_HWLOC:BOOL=OFF
                                      #-DGINKGO_BUILD_OMP:BOOL=OFF
                                      -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}
                                      -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CMAKE_CXX_FLAGS}
                                      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                                      -DCMAKE_C_FLAGS:STRING=${CMAKE_CMAKE_C_FLAGS}
                                      -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
                                      -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_CMAKE_Fortran_FLAGS}
                                      -DCMAKE_INSTALL_PREFIX:PATH=${Ginkgo_install_dir}
                                      -DCMAKE_INSTALL_RPATH:PATH=${Ginkgo_install_dir}/lib
                                      -DCMAKE_INSTALL_NAME_DIR:PATH=${Ginkgo_install_dir}/lib
                                      -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}

                    # -- Build
                    BINARY_DIR       ${Ginkgo_build_dir}        # Build directory 
                    BUILD_COMMAND    $(MAKE)                      # $(MAKE) enables parallel builds through make
                    BUILD_IN_SOURCE  ${Ginkgo_BUILD_IN_SOURCE}  # Flag for in source builds
                    # -- Install
                    INSTALL_DIR      ${Ginkgo_install_dir}      # Install directory
                    # -- Output control
                    ${Ginkgo_logging_args}
            )

# --- Useful variables for packages that depends on Ginkgo
global_set(Ginkgo_INSTALL_PREFIX "${Ginkgo_install_dir}")