# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Process Kernel
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
# include(PrintVariable)
# include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   CHEMPK_SOURCE_DIR
#   CHEMPK_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
#project(CHEMPK)

# Chemistry stuff.
get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)
get_property(GEOCHEM_LIBS GLOBAL PROPERTY GEOCHEM_LIBS)
include_directories(${CHEM_INCLUDES_DIR})

# External (TPL) include directories
# Although it does not call HDF5 directly 
# it does need HDF5 for the Vis.hpp file 
#include_directories(${Epetra_INCLUDE_DIRS})
#include_directories(${HDF5_INCLUDE_DIRS})
# Amanzi include directories
#include_directories(${DBC_SOURCE_DIR})
#include_directories(${ATK_SOURCE_DIR})
#include_directories(${PKS_SOURCE_DIR})
#include_directories(${GEOCHEM_SOURCE_DIR})
#include_directories(${GEOMETRY_SOURCE_DIR})
#include_directories(${MESH_SOURCE_DIR})
#include_directories(${MFUNCS_SOURCE_DIR})
#include_directories(${MESH_DATA_SOURCE_DIR})
#include_directories(${DATA_STRUCTURES_SOURCE_DIR})
#include_directories(${FUNCS_SOURCE_DIR})
#include_directories(${STATE_SOURCE_DIR})
#include_directories(${OUTPUT_SOURCE_DIR})

include_directories(${Amanzi_TPL_MSTK_INCLUDE_DIRS})

#
# Chemistry registrations
#
#register_evaluator_with_factory(
#  HEADERFILE Amanzi_PK_reg.hh
#  LISTNAME   PKS_CHEMISTRY_REGISTRATIONS )

#if (ENABLE_ALQUIMIA)
#register_evaluator_with_factory(
#  HEADERFILE Alquimia_PK_reg.hh
#  LISTNAME   PKS_CHEMISTRY_REGISTRATIONS )
#endif()

#generate_evaluators_registration_header( 
#  HEADERFILE pks_chemistry_registration.hh   
#  LISTNAME   PKS_CHEMISTRY_REGISTRATIONS
#  INSTALL    True )

#
# Install Targets
#
set(chem_inc_files 
#                   Alquimia_PK_ATS.hh
                   Amanzi_PK_ATS.hh
                   Chemistry_PK_ATS.hh
)

#
# Library: chemistry_pk
#
set(chem_src_files) # Chemistry_PK.cc)

list(APPEND chem_src_files Chemistry_PK.cc Amanzi_PK.cc)
if (ENABLE_ALQUIMIA)
   list(APPEND chem_src_files Alquimia_PK.cc)
endif()

# #add_install_include_file(${chem_inc_files})

# add_library(pk_chemistry
#                     SOURCE ${chem_src_files}
#                     HEADERS ${chem_inc_files}
#                     LINK_LIBS ${GEOCHEM_LIBS} pks ${Epetra_LIBRARIES})

 add_library(pk_chemistry
             ${chem_src_files}
)

# install(TARGETS pk_chemistry DESTINATION lib)