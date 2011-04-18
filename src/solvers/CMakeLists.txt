# -*- mode: cmake -*-

#
#  Amanzi
#   Solvers
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(AddParallelTest)

#
# Define a project name
# After this command the following varaibles are defined
#   SOLVERS_SOURCE_DIR
#   SOLVERS_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(SOLVERS)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${NOX_INCLUDE_DIRS})

# Need to move this to roo level -- lpritch
add_definitions("-DUSE_MPI")

#
# Library: solvers
#
file(GLOB solvers_source_files "*.cpp")
add_library(solvers ${solvers_source_files})
set(solvers_tpl_libs ${DBC_LIBRARIES} ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${NOX_LIBRARIES})
list(REMOVE_DUPLICATES solvers_tpl_libs)
target_link_libraries(solvers error_handling ${solvers_tpl_libs})

if (BUILD_TESTS)

endif()
