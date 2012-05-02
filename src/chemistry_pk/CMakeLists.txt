# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Process Kernel
#

#
# Define a project name
# After this command the following varaibles are defined
#   CHEMPK_SOURCE_DIR
#   CHEMPK_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(CHEMPK)

get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)
get_property(GEOCHEM_LIBS GLOBAL PROPERTY GEOCHEM_LIBS)

# Chemistry include
include_directories(${CHEM_INCLUDES_DIR})

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})


# External (TPL) include directories
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${ASCEMIO_INCLUDE_DIR})


add_library(chemistry_pk chemistry_state.cc chemistry_pk.cc)
target_link_libraries(chemistry_pk ${GEOCHEM_LIBS} ${Epetra_LIBRARIES})
add_install_library(chemistry_pk)

if (BUILD_TESTS)
  add_subdirectory(tests_integration)
endif()
