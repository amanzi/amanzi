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
if ( CCSE_FOUND )
  include_directories(${CCSE_INCLUDE_DIRS})
endif()  

include_directories(${Boost_INCLUDE_DIRS})

#
# interface library
#
file(GLOB _H_inc_files "*.H")
file(GLOB _hh_inc_files "*.hh")
file(GLOB _h_inc_files "*.h")
set(interface_inc_files ${_H_inc} ${_hh_inc_files} ${_h_inc_files})
add_amanzi_library(interface
                      SOURCE 
                           InputParser.H InputParser.cpp
                           InputParserIS.hh InputParserIS.cc
                           ObservationData.H Simulator.H
                      HEADERS ${interface_inc_files}
                      LINK_LIBS error_handling geometry ${Teuchos_LIBRARIES})
