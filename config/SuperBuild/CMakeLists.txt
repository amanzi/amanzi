#  -*- mode: cmake -*-

#
# Amanzi SuperBuild
#       Amanzi TPL Build Targets
#
project(SuperBuild)

# Build the docs
add_subdirectory(doc)

# update the module path
set(CMAKE_MODULE_PATH ${SuperBuild_SOURCE_DIR} ${CMAKE_MODULE_PATH})
include(Build_CURL) 
include(Build_zlib)
