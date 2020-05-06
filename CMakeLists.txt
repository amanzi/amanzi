# -*- mode: cmake -*-

#
# Amanzi 
#       Root CMakeLists.txt file
#

# Require cmake version same as the version used to buildg TPLs 
#(see SuperBuild/CMakeLists.txt)
cmake_minimum_required(VERSION 3.6.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable testing, ctest needs this
# all add_test commands are ignored unless this is called!
enable_testing()

# Define the project name
# This command will define
#     AMANZI_SOURCE_DIR
#     AMANZI_BINARY_DIR
project(AMANZI)

message(STATUS "IN TOP CMAKELIST --------------------------------------##################")

message(STATUS "\n-------------------------------------------")
message(STATUS "\n-- CMake: Configuring Amanzi build/install.\n--")
message(STATUS "----------------------------------------")

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


if (SPACKAGE_BUILD)
  set(MSTK_VERSION ON)
  set(MSTK_VERSION_MAJOR 3) 
  set(MSTK_VERSION_MINOR 0)
  # TODO generate the version file 
  set(TPL_VERSIONS_INCLUDE_FILE ${CMAKE_BINARY_DIR}/tpl_versions.h)
  file(WRITE ${TPL_VERSIONS_INCLUDE_FILE} "")
  install(FILES ${CMAKE_BINARY_DIR}/tpl_versions.h DESTINATION include)
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_MPI")
  #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_MPI")
endif()

# Code version defined
include(AmanziVersion)

# Addtional build options
include(AmanziOptions)

message("ETC: Trilinos = ${Trilinos_INSTALL_PREFIX}")


# Find required Amanzi TPL
include(AmanziTPL)

# shared libraries tweaks: enforcing absolute path
if (BUILD_SHARED_LIBS)
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_SKIP_INSTALL_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE) 
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_RPATH}")
else()
  SET(CMAKE_SKIP_INSTALL_RPATH ON CACHE BOOL "Turn off for static install." FORCE)
  SET(CMAKE_SKIP_RPATH ON CACHE BOOL "Turn off for static install." FORCE)
endif()

#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

# Change CMAKE COMPILER if KOKKOS used
if(${ENABLE_CUDA})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -Wno-deprecated-declarations -lineinfo \
    -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored \
    -Xcudafe --diag_suppress=conversion_function_not_usable \
    -Xcudafe --diag_suppress=cc_clobber_ignored \
    -Xcudafe --diag_suppress=code_is_unreachable \
    -Xcudafe --diag_suppress=declared_but_not_referenced \
    -Xcudafe --diag_suppress=integer_sign_change \
    -Xcudafe --diag_suppress=partial_override \
    -Xcudafe --diag_suppress=virtual_function_decl_hidden \
    -Xcudafe --diag_suppress=probable_guiding_friend \
    -Xcudafe --display_error_number \
    --expt-extended-lambda")
  message(STATUS "NVCC_WRAPPER_DEFAULT_COMPILER: $ENV{NVCC_WRAPPER_DEFAULT_COMPILER}")
  message(STATUS "CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER}")
endif()

# Check the mesh framework choice
if ((ENABLE_Unstructured)  AND
    (NOT ENABLE_STK_Mesh)  AND
    (NOT ENABLE_MOAB_Mesh) AND
    (NOT ENABLE_MSTK_Mesh))
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
  find_package(PythonInterp)
  add_subdirectory(test_suites) 
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
    if (NOT UNIX)
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
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
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
