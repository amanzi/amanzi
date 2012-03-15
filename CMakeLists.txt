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
set(CMAKE_MODULE_PATH 
    ${AMANZI_MODULE_PATH}
    ${AMANZI_MODULE_PATH}/Modules
    ${AMANZI_MODULE_PATH}/Utils)

# Code version defined
include(AmanziVersion)

# Addtional build options
include(AmanziOptions)

# Add the TPL build targets
option(ENABLE_SuperBuild "Create TPL build targets" OFF)
if ( ENABLE_SuperBuild )
  add_subdirectory(config/SuperBuild)
endif()

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


# A property for accumulating the a global amanzi link line.
set_property(GLOBAL PROPERTY AMANZI_LINK_LINE "-L${CMAKE_INSTALL_PREFIX}/lib")

# A property for accumulating amanzi library targets
set_property(GLOBAL PROPERTY AMANZI_LIBRARY_TARGETS)
set(AMANZI_LINK_LINE_FILE "${AMANZI_BINARY_DIR}/link_line")  # A filename to write link-line to.
include(InstallManager)

# Source files for all binaries and libraries found under src
add_subdirectory(src)

# Define a list of all enabled TPLs. Must do this AFTER all the 
# CMakelist.txt files have been processed!
get_property(AMANZI_ENABLED_TPLS GLOBAL PROPERTY PACKAGES_FOUND)

# Make the exports only _after_ doing the build
create_exports()

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

