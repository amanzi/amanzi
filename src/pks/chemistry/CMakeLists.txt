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
include_directories(${PKS_SOURCE_DIR})
include_directories(${GEOCHEM_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
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

register_evaluator_with_factory(
  HEADERFILE Amanzi_PK_reg.hh
  LISTNAME   PKS_CHEMISTRY_REGISTRATIONS )

register_evaluator_with_factory(
  HEADERFILE Alquimia_PK_reg.hh
  LISTNAME   PKS_CHEMISTRY_REGISTRATIONS )

generate_evaluators_registration_header( 
  HEADERFILE pks_chemistry_registration.hh   
  LISTNAME   PKS_CHEMISTRY_REGISTRATIONS
  INSTALL    True )


set(chem_inc_files 
                   Alquimia_PK.hh
                   Alquimia_PK_reg.hh
                   Amanzi_PK.hh
                   Amanzi_PK_reg.hh
                   Chemistry_PK.hh
                   pks_chemistry_registration.hh
)

set(chem_src_files)
if (ENABLE_Unstructured)
  list(APPEND chem_src_files Chemistry_PK.cc Amanzi_PK.cc)
  if (ENABLE_ALQUIMIA)
    list(APPEND chem_src_files Alquimia_PK.cc)
  endif()
endif()

add_amanzi_library(chemistry_pk
                   SOURCE ${chem_src_files}
                   HEADERS ${chem_inc_files}
                   LINK_LIBS ${GEOCHEM_LIBS} ${Epetra_LIBRARIES})

if (BUILD_TESTS AND ENABLE_Unstructured)
  add_subdirectory(tests_integration)
endif()
