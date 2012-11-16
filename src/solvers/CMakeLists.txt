# -*- mode: cmake -*-

#
#  Amanzi
#   Solvers
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

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
include_directories(${DBG_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${NOX_INCLUDE_DIRS})

#
# Library: solvers
#
file(GLOB solvers_source_files "*.cpp")
file(GLOB solver_inc_files "*.H")
set(solvers_tpl_libs ${DBC_LIBRARIES} ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${NOX_LIBRARIES})
list(REVERSE solvers_tpl_libs)
list(REMOVE_DUPLICATES solvers_tpl_libs)
list(REVERSE solvers_tpl_libs)
add_amanzi_library(solvers
                   SOURCE ${solvers_source_files} HEADERS ${solvers_inc_files}
                   LINK_LIBS error_handling ${solvers_tpl_libs})

if (BUILD_TESTS)

endif()
