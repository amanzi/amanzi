# -*- mode: cmake -*-

#
#  Amanzi
#    Error Handling A.K.A Design by contract (DBC), Throw, Catch, etc.
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)
include(InstallManager)

#
# Define a project name
# After this command the following varaibles are defined
#   DBC_SOURCE_DIR
#   DBC_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(DBC)

#
# Define a the list of include files
#  will use this if/when we define an install make target
file(GLOB dbc_inc_files "*.hh")

#
# Define source files
file(GLOB dbc_lib_files "*.cc")

# Define the library 
# add_library(<libaname> file1 file2 ... )
add_library(error_handling ${dbc_lib_files})
link_list_add(error_handling)

#
# Install Targets
#
add_install_include_file(${dbc_inc_files})
add_install_library(error_handling)

if (BUILD_TESTS)
   
  # Get the UnitTest libraries
  if ( NOT UnitTest_FOUND ) 
      find_package(UnitTest REQUIRED)
  endif()  

  # Add UnitTest include directory
  include_directories(${UnitTest_INCLUDE_DIRS})

  # Define the test executable
  # File containing the 'main' symbol should be listed first
  add_executable(test_error_handling 
                 test/Main.cc test/test_exceptions.cc test/test_dbc.cc test/test_errors.cc)
  
  # Define the links this binary needs
  target_link_libraries(test_error_handling error_handling ${UnitTest_LIBRARIES})

  # Add the executable to the 'make test' target
  # this command does nothing if enable_testing() is not 
  # called in a parent CMakeLists.txt file 
  add_amanzi_test(error test_error_handling KIND unit)

endif()



