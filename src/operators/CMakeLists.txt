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

include_directories(${ATS_SOURCE_DIR}/src/operators)
add_subdirectory(advection)

# elliptic operators?