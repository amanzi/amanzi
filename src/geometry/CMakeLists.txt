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

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})

#
# geometry library
#
add_library(geometry Geometry.cc Point.hh 
                     Region.hh Region.cc
                     RectangularRegion.hh RectangularRegion.cc)
target_link_libraries(geometry error_handling)

#
# Install Targets
#
add_install_library(geometry)
file(GLOB _inc_files "*.h*")
add_install_include_file(${_inc_files})


