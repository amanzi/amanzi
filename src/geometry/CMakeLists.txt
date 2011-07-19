# -*- mode: cmake -*-

#
#  Amanzi
#    Mesh 
#
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)

#
# Define a project name
# After this command the following varaibles are defined
#   GEOMETRY_SOURCE_DIR
#   GEOMETRY_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(GEOMETRY)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})

#
# geometry library
#
add_library(geometry Geometry.cc Point.hh)
add_install_library(geometry)
link_list_add(geometry)

#
# Install targets
#
add_install_library(geometry)

