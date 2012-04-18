# -*- mode: cmake -*-


#
#  Amanzi Tool Kit (AKT)
#    General utilities
#    This directory has become a dumping ground for things that do not
#    have a well-defined home. Will go away or get renamed. For now we
#    need to keep it around.
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)

#
# Define a project name
# After this command the following varaibles are defined
#   ATK_SOURCE_DIR
#   ATK_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(ATK)

set(atk_inc_files Utils.hh
                  Timer.hh
)

add_library(atk Timer.cc
)

add_install_library(atk)

#
# Install Targets
#
add_install_include_file(${atk_inc_files})













