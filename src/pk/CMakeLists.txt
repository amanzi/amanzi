# -*- mode: cmake -*-

#
#  Amanzi
#    PK class
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   PK_SOURCE_DIR
#   PK_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(PK)

# Need to move this up to the root -- lpritch
add_definitions("-DUSE_MPI")

#
# Install Target
#
file(GLOB PK_inc_files "*.hh")
add_install_include_file(${PK_inc_files})
