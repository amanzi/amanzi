# -*- mode: cmake -*-

#
#  Amanzi
#    Common geochemistry machinery
#

#
# Define a project name
# After this command the following variables are defined
#   GEOCHEM_SOURCE_DIR
#   GEOCHEM_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(GEOCHEM)

# We use Alquimia for chemistry instead of Amanzi's built-in thing, 
# so we augment whatever existing GEOCHEM_LIBS property exists.
get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)
set(CHEM_INCLUDES_DIR ${CHEM_INCLUDES_DIR} ${GEOCHEM_SOURCE_DIR} ${Alquimia_INCLUDE_DIRS})
set_property(GLOBAL PROPERTY CHEM_INCLUDES_DIR "${CHEM_INCLUDES_DIR}")
get_property(GEOCHEM_LIBS GLOBAL PROPERTY GEOCHEM_LIBS)
set(GEOCHEM_SUPPORT_LIBS ${GEOCHEM_LIBS} ${Alquimia_LIBRARIES} ${PFLOTRAN_LIBRARIES} ${CrunchTope_LIBRARIES} ${PETSc_LIBRARIES})
set(GEOCHEM_LIBS ${GEOCHEM_LIBS} geochemistry ${GEOCHEM_SUPPORT_LIBS})
set_property(GLOBAL PROPERTY GEOCHEM_LIBS "${GEOCHEM_LIBS}")

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${CHEM_INCLUDES_DIR})

# External (TPL) include directories
include_directories(${Epetra_INCLUDE_DIRS})

file(GLOB chem_inc_files "*hh")
set(chem_src_files ChemistryEngine.cc)
add_amanzi_library(geochemistry
                   SOURCE ${chem_src_files}
                   HEADERS ${chem_inc_files}
                   LINK_LIBS ${GEOCHEM_SUPPORT_LIBS} ${Epetra_LIBRARIES})

