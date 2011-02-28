# -*- mode: cmake -*-

#
# Amanzi 
#       Root CMakeLists.txt file
#

# Require cmake 2.8 or higher
cmake_minimum_required(VERSION 2.8)

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

# Addtional build options
include(AmanziOptions)

# Find required Amanzi TPL
include(AmanziTPL)

# Source files for all binaries and libraries found under src
add_subdirectory(src)

option(ENABLE_Config_Report "Print out detailed information at the end of a configuration" ON)
if (ENABLE_Config_Report)
    include(AmanziConfigReport)
endif()    
