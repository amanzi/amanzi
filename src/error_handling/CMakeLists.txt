# -*- mode: cmake -*-

#
#  Amanzi
#    Error Handling A.K.A Design by contract (DBC), Throw, Catch, etc.
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)
include(LibraryManager)

#
# Define a project name
# After this command the following varaibles are defined
#   DBC_SOURCE_DIR
#   DBC_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(DBC)

if(AMANZI_ARCH STREQUAL "Summit")
  message(STATUS "SUMMIT ARCH")
  message(STATUS "INCLUDE: ${Kokkos_INCLUDE_DIRS}")
  message(STATUS "INCLUDE: ${CUDA_INCLUDE_DIRS}")
  include_directories(${Kokkos_INCLUDE_DIRS})
  include_directories(${CUDA_INCLUDE_DIRS})
endif()


#
# Library Target: error_handling
#
file(GLOB dbc_lib_files "*.cc")
file(GLOB dbc_inc_files "*.hh")
add_amanzi_library(error_handling SOURCE ${dbc_lib_files} HEADERS ${dbc_inc_files} LINK_LIBS ${Kokkos_LIBRARIES} ${CUDA_LIBRARIES})

if (BUILD_TESTS)
  # Get the UnitTest libraries
  if ( NOT UnitTest_FOUND ) 
      find_package(UnitTest REQUIRED)
  endif()  

  # Add UnitTest include directory
  include_directories(${UnitTest_INCLUDE_DIRS})

  # Create the error test suite
  file(GLOB error_test_files RELATIVE ${DBC_SOURCE_DIR} "${DBC_SOURCE_DIR}/test/test*.cc")
  add_amanzi_test(error test_error_handling
                  KIND unit
		  SOURCE test/Main.cc ${error_test_files}
		  LINK_LIBS error_handling ${UnitTest_LIBRARIES})
endif()



