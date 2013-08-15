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

IF (ENABLE_ALQUIMIA)
  # We use Alquimia for chemistry instead of Amanzi's built-in thing, 
  # so we replace whatever existing GEOCHEM_LIBS property exists.

  # The following code, which expresses the libraries for Alquimia and 
  # PFlotran to CMake, belongs somewhere else, but there's no obvious home
  # for it.

  # PFlotran libraries.
  add_library(pflotranchem STATIC IMPORTED)
  set_target_properties(pflotranchem PROPERTIES IMPORTED_LOCATION "${CMAKE_INSTALL_PREFIX}/tpls/lib/libpflotranchem.a")

  # Alquimia libraries.
  add_library(alquimia_c STATIC IMPORTED)
  set_target_properties(alquimia_c PROPERTIES IMPORTED_LOCATION "${CMAKE_INSTALL_PREFIX}/tpls/lib/libalquimia_c.a")
  add_library(alquimia_cutils STATIC IMPORTED)
  set_target_properties(alquimia_cutils PROPERTIES IMPORTED_LOCATION "${CMAKE_INSTALL_PREFIX}/tpls/lib/libalquimia_cutils.a")
  add_library(alquimia_fortran STATIC IMPORTED)
  set_target_properties(alquimia_fortran PROPERTIES IMPORTED_LOCATION "${CMAKE_INSTALL_PREFIX}/tpls/lib/libalquimia_fortran.a")

  set(ALQUIMIA_LIBRARIES alquimia_c;alquimia_cutils;alquimia_fortran;pflotranchem;${PETSC_LIBRARIES})
  set(GEOCHEM_LIBS ${ALQUIMIA_LIBRARIES})
  set_property(GLOBAL PROPERTY GEOCHEM_LIBS "${GEOCHEM_LIBS}")
  set(CHEM_INCLUDES_DIR ${ALQUIMIA_INCLUDE_DIRS})
  set_property(GLOBAL PROPERTY CHEM_INCLUDES_DIR "${CHEM_INCLUDES_DIR}")
ELSE()
  # Just use what's already there.
  get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)
  get_property(GEOCHEM_LIBS GLOBAL PROPERTY GEOCHEM_LIBS)
ENDIF()

# Chemistry include
include_directories(${CHEM_INCLUDES_DIR})

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
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
set(chem_src_files Chemistry_State.cc)
if (ENABLE_ALQUIMIA)
  list(APPEND chem_src_files alquimia_chemistry_pk.cc)
else()
  list(APPEND chem_src_files chemistry_pk.cc)
endif()
add_amanzi_library(chemistry_pk
                   SOURCE ${chem_src_files}
                   HEADERS ${chem_inc_files}
                   LINK_LIBS ${GEOCHEM_LIBS} ${Epetra_LIBRARIES})

if (BUILD_TESTS)
  if (NOT ENABLE_ALQUIMIA)
    add_subdirectory(tests_integration)
  endif()
endif()
