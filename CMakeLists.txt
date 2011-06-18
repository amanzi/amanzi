# -*- mode: cmake -*-

#
# Amanzi 
#       Root CMakeLists.txt file
#

# Require cmake 2.8 or higher
cmake_minimum_required(VERSION 2.8.2)

if ("${CMAKE_PATCH_VERSION}" GREATER 3)
  set(ADJUST_POLICY "True")
endif()

if (${ADJUST_POLICY})
  cmake_policy(SET CMP0017 NEW)
endif()

# Enable testing, ctest needs this
# all add_test commands are ignored unless this is called!
enable_testing()

# Define the project name
# This command will define
#     AMANZI_SOURCE_DIR
#     AMANZI_BINARY_DIR
project(AMANZI)

# Useful variables pointing to directories in the source tree
set(AMANZI_SOURCE_SRC_DIR      "${AMANZI_SOURCE_DIR}/src")
set(AMANZI_SOURCE_TOOLS_DIR    "${AMANZI_SOURCE_DIR}/tools")
set(AMANZI_SOURCE_EXAMPLES_DIR "${AMANZI_SOURCE_DIR}/examples")
set(AMANZI_PYTHON_DIR          "${AMANZI_SOURCE_DIR}/tools/py_lib")

# Set the module search path so find_package and include commands
# can locate files in <root source tree>/tools/cmake
set(AMANZI_MODULE_PATH "${AMANZI_SOURCE_TOOLS_DIR}/cmake")
set(CMAKE_MODULE_PATH  ${AMANZI_MODULE_PATH})

# Code version defined
include(AmanziVersion)

#
# Install locations defined
#
# Use CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT
# to override the /usr/local as the default prefix.
# This is not working! 
set(DEFAULT_INSTALL_PREFIX "$ENV{HOME}/amanzi-${AMANZI_VERSION}")
if (CMAKE_INSTALL_PREFIX_INTIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${DEFAULT_INSTALL_PREFIX}" 
                         CACHE STRING "Install path prefix" FORCE)
    message(STATUS "Using default install prefix")
else()
    message(STATUS "User defined install prefix  "
	           "${CMAKE_INSTALL_PREFIX}")
endif()    
set(AMANZI_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
set(AMANZI_INSTALL_INCLUDE_PATH "${AMANZI_INSTALL_PREFIX}/include")
set(AMANZI_INSTALL_LIB_PATH     "${AMANZI_INSTALL_PREFIX}/lib")
set(AMANZI_INSTALL_BIN_PATH     "${AMANZI_INSTALL_PREFIX}/bin")
message(STATUS "Amanzi install prefix ${AMANZI_INSTALL_PREFIX}")
message(STATUS "\t include files ${AMANZI_INSTALL_INCLUDE_PATH}")
message(STATUS "\t libraries ${AMANZI_INSTALL_LIB_PATH}")
message(STATUS "\t binaries ${AMANZI_INSTALL_BIN_PATH}")




# Addtional build options
include(AmanziOptions)

# Find required Amanzi TPL
include(AmanziTPL)

# Check the mesh framework choice
if ( (NOT ENABLE_STK_Mesh)  AND 
     (NOT ENABLE_MOAB_Mesh) AND
     (NOT ENABLE_MSTK_Mesh)
   )
    message(FATAL_ERROR "Missing a mesh framework\n"
                        "Please enable at least one of the following mesh frameworks\n"
                        "-D ENABLE_MOAB_Mesh:BOOL=ON\n"
                        "-D ENABLE_STK_Mesh:BOOL=ON\n"
                        "-D ENABLE_MSTK_Mesh:BOOL=ON\n")
endif()    

# Source files for all binaries and libraries found under src
add_subdirectory(src)

# Create the export Makefile
include(InstallManager)
create_export_makefile()

# Write out the configuration log file and dump to the screen if requested
option(ENABLE_Config_Report "Print out detailed information at the end of a configuration")
set(AMANZI_CONFIG_LOG "${AMANZI_BINARY_DIR}/amanzi-config.log"
    CACHE string "Amanzi configuration log file")
include(AmanziConfigReport)
if (ENABLE_Config_Report)
    set(cat_exec "cat")
    if (WIN32)
        if( NOT UNIX)
            set(cat_exec "type")
        endif(NOT UNIX)
    endif(WIN32)

    execute_process(COMMAND "${cat_exec}" "${AMANZI_CONFIG_LOG}" OUTPUT_VARIABLE config_output)
    print_variable(cat_exec)
    print_variable(AMANZI_CONFIG_LOG)
    print_variable(config_output)
    message(STATUS "********************************************************************************")
    message(STATUS "begin configuration output --\n${config_output}")
    message(STATUS "end configuration output --")
    message(STATUS "********************************************************************************")

endif()


