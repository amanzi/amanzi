# -*- mode: cmake -*-

#
#  ATS
#   Operators
#

#
# Define a project name
# After this command the following varaibles are defined
#   OPERATOR_SOURCE_DIR
#   OPERATOR_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.

add_subdirectory(advection)
add_subdirectory(divgrad)

# elliptic operators?