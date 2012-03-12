#  -*- mode: cmake -*-

#
# Build TPL: ExodusII 
#    
define_external_project(ExodusII 
                        TARGET exodusii
                        DEPENDS netcdf)

# ########################################################################### #
# Build the configure command
# ########################################################################### #

# Need to define variables used in the configure script
include(BuildWhitespaceString)
build_whitespace_string(common_cmake_args
                       ${Amanzi_CMAKE_COMPILER_ARGS}
		       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})

# Build the configure script
set(ExodusII_sh_configure ${ExodusII_prefix_dir}/exodusii-configure-step.sh)
configure_file(${SuperBuild_BUILD_FILES_DIR}/exodusii-configure-step.sh.in
               ${ExodusII_sh_configure}
        @ONLY)

# Configure the CMake command file
set(ExodusII_cmake_configure ${ExodusII_prefix_dir}/exodusii-configure-step.cmake)
configure_file(${SuperBuild_BUILD_FILES_DIR}/exodusii-configure-step.cmake.in
               ${ExodusII_cmake_configure}
        @ONLY)
set(ExodusII_CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${ExodusII_cmake_configure})  

# ########################################################################### #
# Add External Project
# ########################################################################### #
ExternalProject_Add(${ExodusII_target}
    DEPENDS netcdf
    ${ExodusII_ep_directory_args}
    ${ExodusII_url_args}
    CONFIGURE_COMMAND ${ExodusII_CONFIGURE_COMMAND}
    ${ExodusII_logging_opts}
)

# ########################################################################### #
# Build nemesis
# ########################################################################### #

# Configure the build script
set(NEMESIS_sh_build ${ExodusII_prefix_dir}/nemesis-build.sh)
configure_file(${SuperBuild_BUILD_FILES_DIR}/nemesis-build.sh.in
               ${NEMESIS_sh_build}
        @ONLY)

# Configure the CMake command file
set(NEMESIS_cmake_build ${ExodusII_prefix_dir}/nemesis-build.cmake)
configure_file(${SuperBuild_BUILD_FILES_DIR}/nemesis-build.cmake.in
               ${NEMESIS_cmake_build}
            @ONLY)
set(NEMESIS_BUILD_COMMAND ${CMAKE_COMMAND} -P ${NEMESIS_cmake_build}) 

ExternalProject_Add_Step(${ExodusII_target} nemesis
                         COMMAND ${NEMESIS_BUILD_COMMAND}
                         COMMENT "Building nemesis (ExodusII extension)"
                         DEPENDEES install
                         WORKING_DIRECTORY ${ExodusII_prefix}
                         LOG TRUE)
