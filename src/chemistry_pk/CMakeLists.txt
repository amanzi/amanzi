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

# Chemistry stuff.
get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)
get_property(GEOCHEM_LIBS GLOBAL PROPERTY GEOCHEM_LIBS)
include_directories(${CHEM_INCLUDES_DIR})

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${GEOCHEM_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})


# External (TPL) include directories
# Although it does not call HDF5 directly 
# it does need HDF5 for the Vis.hpp file 
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})

file(GLOB chem_inc_files "*hh")
set(chem_src_files Chemistry_State.cc;chemistry_pk.cc)
if (ENABLE_ALQUIMIA)
  list(APPEND chem_src_files alquimia_chemistry_pk.cc)
else()
endif()
add_amanzi_library(chemistry_pk
                   SOURCE ${chem_src_files}
                   HEADERS ${chem_inc_files}
                   LINK_LIBS ${GEOCHEM_LIBS} ${Epetra_LIBRARIES})

if (BUILD_TESTS)
  add_subdirectory(tests_integration)
endif()
