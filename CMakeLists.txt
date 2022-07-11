# -*- mode: cmake -*-

#
# Amanzi 
#       Root CMakeLists.txt file
#

# Require cmake version same as the version used to buildg TPLs 
#(see SuperBuild/CMakeLists.txt)
cmake_minimum_required(VERSION 3.13.3)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

# Enable testing, ctest needs this
# all add_test commands are ignored unless this is called!
enable_testing()

# Define the project name
# This command will define
#     AMANZI_SOURCE_DIR
#     AMANZI_BINARY_DIR
project(AMANZI)

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

# Append user supplied flags to base CMAKE_<LANL>_FLAGS and remove duplicates

string(CONCAT CMAKE_C_FLAGS_TMP ${TPL_CMAKE_C_FLAGS} " " ${CMAKE_C_FLAGS})
string(CONCAT CMAKE_CXX_FLAGS_TMP ${TPL_CMAKE_CXX_FLAGS} " " ${CMAKE_CXX_FLAGS})
string(CONCAT CMAKE_Fortran_FLAGS_TMP ${TPL_CMAKE_Fortran_FLAGS} " " ${CMAKE_Fortran_FLAGS})

separate_arguments(CMAKE_C_FLAGS_TMP)
separate_arguments(CMAKE_CXX_FLAGS_TMP)
separate_arguments(CMAKE_Fortran_FLAGS_TMP)

list(REMOVE_DUPLICATES CMAKE_C_FLAGS_TMP)
list(REMOVE_DUPLICATES CMAKE_CXX_FLAGS_TMP)
list(REMOVE_DUPLICATES CMAKE_Fortran_FLAGS_TMP)

string(JOIN " " CMAKE_C_FLAGS ${CMAKE_C_FLAGS_TMP})
string(JOIN " " CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_TMP})
string(JOIN " " CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_TMP})

# Code version defined
include(AmanziVersion)

# Addtional build options
include(AmanziOptions)

# Find required Amanzi TPL
include(AmanziTPL)
# The following line allows us to link third-party libraries not explicitly 
# found by the logic in AmanziTPL.cmake.
#
# ETC: this is wrong, should point to TPLs, but not necessarily installed here, meaning that lots of ld warnings show up.
#link_directories(${CMAKE_INSTALL_PREFIX}/tpls/lib)

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

# Check the mesh framework choice
if ((ENABLE_Unstructured)  AND
    (NOT ENABLE_MESH_MOAB) AND
    (NOT ENABLE_MESH_MSTK))
  message(FATAL_ERROR "Missing a mesh framework\n"
                      "Please enable at least one of the following mesh frameworks\n"
                      "-D ENABLE_MESH_MOAB:BOOL=ON\n"
                      "-D ENABLE_MESH_MSTK:BOOL=ON\n")
endif()    

if(ENABLE_OpenMP)
  find_package(OpenMP REQUIRED)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if(ENABLE_CUDA OR ENABLE_UVM) 
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -arch=sm_70 \
    -lineinfo -Wno-deprecated-declarations \
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
    -lcublas -lcudadevrt -lcuda -lcudart -lcusparse \
    --expt-extended-lambda")  
  message(STATUS "NVCC_WRAPPER_DEFAULT_COMPILER: $ENV{NVCC_WRAPPER_DEFAULT_COMPILER}")  
  message(STATUS "CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER}")  
endif()

# A property for accumulating the a global amanzi link line.
set_property(GLOBAL PROPERTY AMANZI_LINK_LINE "-L${CMAKE_INSTALL_PREFIX}/lib")

# A property for accumulating amanzi library targets
set_property(GLOBAL PROPERTY AMANZI_LIBRARY_TARGETS)
set(AMANZI_LINK_LINE_FILE "${AMANZI_BINARY_DIR}/link_line")  # A filename to write link-line to.
include(InstallManager)

# find python3 for use in regression testing
if (ENABLE_Regression_Tests)
  find_package(PythonInterp 3 REQUIRED VERSION)
endif()

# Source files for all binaries and libraries found under src
add_subdirectory(src)

# Examples added as tests
# add_subdirectory(examples) 

# Adding verification tests
# EIB - uncommenting this will include the verification/benchmarking tests
include(CTest)
if (ENABLE_Regression_Tests)
  if (ENABLE_AmanziPhysicsModule)
    add_subdirectory(test_suites)
  endif()
endif()


# Python modules
add_subdirectory(tools/py_lib)

# Define a list of all enabled TPLs. Must do this AFTER all the 
# CMakelist.txt files have been processed!
get_property(AMANZI_ENABLED_TPLS GLOBAL PROPERTY PACKAGES_FOUND)


# Show TPL and current CMake_<LANG>_FLAGS
#
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UC)

message(STATUS "------------------------------------------------")
message(STATUS "CMake Compiler Configuration for TPLs ... ")
message(STATUS "TPL_CMAKE_C_FLAGS       = ${TPL_CMAKE_C_FLAGS}")
message(STATUS "TPL_CMAKE_CXX_FLAGS     = ${TPL_CMAKE_CXX_FLAGS}")
message(STATUS "TPL_CMAKE_Fortran_FLAGS = ${TPL_CMAKE_Fortran_FLAGS}")
message(STATUS "------------------------------------------------")
message(STATUS "")
message(STATUS "CMake Compiler Configuration ... ")
message(STATUS "")
message(STATUS "CMAKE_C_FLAGS       = ${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UC} = ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UC}}")
message(STATUS "CMAKE_CXX_FLAGS     = ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UC} = ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UC}}")
message(STATUS "CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")
message(STATUS "CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UC} = ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UC}}")
message(STATUS "CMAKE_BUILD_TYPE    = ${CMAKE_BUILD_TYPE}")
message(STATUS "------------------------------------------------")

# Make the exports only _after_ doing the build
create_exports()

option(ENABLE_Config_Report "Print out detailed information at the end of a configuration")
set(AMANZI_CONFIG_LOG "${AMANZI_BINARY_DIR}/amanzi-config.log"
    CACHE STRING "Amanzi configuration log file")
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
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYRIGHT.md")
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
