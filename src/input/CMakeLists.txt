# -*- mode: cmake -*-

#
#  Amanzi
#    Input 
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(AddParallelTest)

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



# Add examples
option(ENABLE_INPUT_EXAMPLES "Build input examples" OFF)
if (ENABLE_INPUT_EXAMPLES)
    add_subdirectory(examples)
endif()    

