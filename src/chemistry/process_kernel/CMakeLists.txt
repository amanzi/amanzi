# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Process Kernel
#

# Chemistry include
include_directories(${CHEM_INCLUDES_DIR})

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_BASE_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Epetra_INCLUDE_DIRS})

set(GEOCHEM_LIBS
	geochemrxns
	geochembase
	geochemsolvers
	geochemutil
)

add_library(chemistry Chemistry_PK.cpp Chemistry_State.cpp)
target_link_libraries(chemistry ${GEOCHEM_LIBS} ${Epetra_LIBRARIES})

