# -*- mode: cmake -*-

#
# Amanzi 
#       Root CMakeLists.txt file
#

# Require cmake 2.8.8 or higher
cmake_minimum_required(VERSION 2.8.8)
set(CMAKE_CTEST_COMMAND ${CMAKE_CTEST_COMMAND})  

if ("${CMAKE_PATCH_VERSION}" GREATER 11)
  cmake_policy(SET CMP0022 NEW)
endif()

if ("${CMAKE_MAJOR_VERSION}" GREATER 2)
  cmake_policy(SET CMP0037 OLD)
  cmake_policy(SET CMP0048 OLD)
endif()

if (CMAKE_VERSION VERSION_LESS "3.1")
  set (CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++ ${CMAKE_CXX_FLAGS}")
else ()
  set(CMAKE_CXX_STANDARD 11)
endif ()

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

# Find required Amanzi TPL
include(AmanziTPL)
# The following line allows us to link third-party libraries not explicitly 
# found by the logic in AmanziTPL.cmake.
link_directories(${CMAKE_INSTALL_PREFIX}/tpls/lib)

# Check the mesh framework choice
if ( (ENABLE_Unstructured)  AND
     (NOT ENABLE_STK_Mesh)  AND
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

# Examples added as tests
# add_subdirectory(examples) 

# Adding verification tests
# EIB - uncommenting this will include the verification/benchmarking tests
include(CTest)
if (ENABLE_Regression_Tests)
    add_subdirectory(testing) 
endif()


# Python modules
add_subdirectory(tools/py_lib)

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

#
#  Create source package
# 
#  - use the TGZ generator
#  - use our existing CMake / Mercurial hooks to get version information
#  - include the generated amanzi_version.hh in the package
#
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Amanzi flow and reactive transport simulator.")
SET(CPACK_PACKAGE_VENDOR "Amanzi Development Team (LANL, LBNL, PNNL)")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYRIGHT")
SET(CPACK_PACKAGE_VERSION_MAJOR ${AMANZI_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${AMANZI_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${AMANZI_VERSION_PATCH}_${AMANZI_VERSION_HASH})
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "amanzi-${AMANZI_VERSION}")
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_IGNORE_FILES ".hg;.hgtags")
SET(CPACK_SOURCE_INSTALLED_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR};/") 
# Copy over extra files: stash them first and then they are copied as part of "make package_source"
LIST(APPEND CPACK_SOURCE_INSTALLED_DIRECTORIES ${CMAKE_CURRENT_BINARY_DIR}/extras .)

INCLUDE(CPack)
