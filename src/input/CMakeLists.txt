# -*- mode: cmake -*-

#
#  Amanzi
#    Input 
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following variables are defined
#   XXX_SOURCE_DIR
#   XXX_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(INPUT)


# Define include directories to build any binary or library
# in this directory. 

# Amanzi include directories

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})


#
# Library: input
#
set(input_src_files input_base_class.cc simple_mesh_input.cc
                    transport_input.cc)
add_library(input ${input_src_files})
target_link_libraries(input ${Teuchos_LIBRARIES})

#
# Library: input_test
#
# This is test input library to exercise language interfaces
# 
set(test_input_src test_input.cc)
add_library(input_test ${test_input_src})
target_link_libraries(input_test ${Teuchos_LIBRARIES})

# Add examples
option(ENABLE_INPUT_EXAMPLES "Build input examples" OFF)
if (ENABLE_INPUT_EXAMPLES)
    add_subdirectory(examples)
endif()    

