# -*- mode: cmake -*-

#
#  Amanzi
#    Interface Platform
#
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)

#
# Define a project name
# After this command the following varaibles are defined
#   PLATFORM_INTERFACE_SOURCE_DIR
#   PLATFORM_INTERFACE_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(PLATFORM_INTERFACE)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${CCSE_INCLUDE_DIRS})

#
# interface library
#
add_library(interface InputParser.H InputParser.cpp
                      InputParserIS.hh InputParserIS.cc
                      ObservationData.H Simulator.H)
target_link_libraries(interface error_handling geometry ${Teuchos_LIBRARIES}) 		    
add_install_library(interface)

#
# Add include files to the installation target
#
file(GLOB _inc_files "*.H *.hh")
add_install_include_file(${_inc_files})
