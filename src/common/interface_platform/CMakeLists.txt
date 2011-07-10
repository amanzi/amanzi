# -*- mode: cmake -*-
#
# Define a project name
# After this command the following varaibles are defined
#   PLATFORM_INTERFACE_SOURCE_DIR
#   PLATFORM_INTERFACE_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(PLATFORM_INTERFACE)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIRS})

if (BUILD_TESTS)
endif()





