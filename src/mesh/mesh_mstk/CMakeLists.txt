# -*- mode: cmake -*-

#
#  Amanzi
#    MESH_MSTK frameworks
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   MESH_MSTK_SOURCE_DIR
#   MESH_MSTK_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MESH_MSTK)

# External (TPL) include directories
include_directories(${MSTK_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${DBC_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})


# Need this define. Errors from MSTK include files 
# about MPI_COMM_WORLD. --lpritch
add_definitions("-DMSTK_HAVE_MPI")

if (WITH_MSTK_2_20rc1_OR_NEWER)
  add_definitions("-DMSTK_2_20rc1_OR_NEWER")
endif ()

if (WITH_MSTK_2_21rc1_OR_NEWER)
  add_definitions("-DMSTK_2_21rc1_OR_NEWER")
endif ()

if (WITH_MSTK_3_00_OR_NEWER)
  add_definitions("-DMSTK_3_00_OR_NEWER")
endif ()

#
# Library: mesh_mstk
#
set(mesh_mstk_files Mesh_MSTK.cc Mesh_MSTK_deform.cc) 
file(GLOB mesh_mstk_inc_files "*.hh")
add_amanzi_library(mesh_mstk
                   SOURCE ${mesh_mstk_files} HEADERS ${mesh_mstk_inc_files}
                   LINK_LIBS geometry mesh atk ${MSTK_LIBRARIES} ${Epetra_LIBRARIES}
		             ${Teuchos_LIBRARIES} ${Zoltan_LIBRARIES} ${HDF5_LIBRARIES})

if (BUILD_TESTS)
    
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Need to copy files for the tests 
    if (NOT (${MESH_MSTK_SOURCE_DIR} EQUAL ${MESH_MSTK_BINARY_DIR}))
        file(GLOB DataFiles "${MESH_MSTK_SOURCE_DIR}/test/*.xml"
                            "${MESH_MSTK_SOURCE_DIR}/test/*.txt"
                            "${MESH_MSTK_SOURCE_DIR}/test/*.exo"
                            "${MESH_MSTK_SOURCE_DIR}/test/*.par.*")
        file(COPY ${DataFiles} DESTINATION ${MESH_MSTK_BINARY_DIR}/test/)
    endif()

    #Test: extracted surface from cell from volumetric face label
    add_amanzi_test(mesh_mstk_surface_face_inheritance mesh_mstk_surface_face_inheritance
                    KIND unit
                    SOURCE test/Main.cc test/mesh_mstk_label_inheritance.cc
                    LINK_LIBS mesh_mstk geometry ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})
endif()

