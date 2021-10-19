#
#  Amanzi
#    Mesh Extracted framework
#


# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)

#
# Define a project name
# After this command the following varaibles are defined
#   MESH_EXTRACTED_SOURCE_DIR
#   MESH_EXTRACTED_BINARY_DIR
project(MESH_EXTRACTED)

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
# include_directories(${MESH_SIMPLE_SOURCE_DIR})

# if (ENABLE_MESH_MOAB)
#   include_directories(${MESH_MSTK_SOURCE_DIR})
# endif()
# if (ENABLE_MESH_MOAB)
#   include_directories(${MESH_MOAB_SOURCE_DIR})
# endif()

# External (TPL) directories
include_directories(${Epetra_INCLUDE_DIRS})
# include_directories(${MOAB_INCLUDE_DIRS})

#
# Library: mesh
#
set(mesh_extracted_inc_files
    BlockMapUtils.hh
    MeshExtractedManifold.hh
)

add_amanzi_library(mesh_extracted
                   SOURCE BlockMapUtils.cc
                          #MeshExtractedManifold.cc
                          MeshFrameworkColumn.cc
                          #MeshSurfaceCell.cc
                   HEADERS ${mesh_extracted_inc_files}
                   LINK_LIBS mesh geometry ${Epetra_LIBRARIES})

if (BUILD_TESTS)
  include_directories(${UnitTest_INCLUDE_DIRS})

  include_directories(${MESH_SIMPLE_SOURCE_DIR})
  list(APPEND mesh_frameworks mesh_simple)

  if (ENABLE_MESH_MOAB)
    include_directories(${MESH_MOAB_SOURCE_DIR})
    include_directories(${MOAB_INCLUDE_DIR})
    list(APPEND mesh_frameworks mesh_moab)
    add_definitions("-DHAVE_MESH_MOAB")
  endif()

  if (ENABLE_MESH_MSTK)
    include_directories(${MESH_MSTK_SOURCE_DIR})
    include_directories(${MSTK_INCLUDE_DIR})
    list(APPEND mesh_frameworks mesh_mstk)
    add_definitions("-DHAVE_MESH_MSTK -DMSTK_HAVE_MPI")
  endif()

  include_directories(${MESH_EXTRACTED_SOURCE_DIR})

  # Need to copy files for the tests 
  if (NOT (${MESH_EXTRACTED_SOURCE_DIR} EQUAL ${MESH_EXTRACTED_BINARY_DIR}))
    file(GLOB DataFiles "${MESH_EXTRACTED_SOURCE_DIR}/test/*.exo"
                           "${MESH_EXTRACTED_SOURCE_DIR}/test/*.xml")
    file(COPY ${DataFiles} DESTINATION ${MESH_EXTRACTED_BINARY_DIR}/test/)
  endif()

  # add_amanzi_test(mesh_extracted_fracture mesh_extracted_fracture
  #                 KIND unit
  #                 SOURCE test/Main.cc test/mesh_extracted_fracture.cc
  #                 LINK_LIBS ${mesh_frameworks} mesh_extracted mesh_audit geometry ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})
  # add_amanzi_test(mesh_extracted_fracture_3pe mesh_extracted_fracture NPROCS 3 KIND unit)

  # add_amanzi_test(mesh_surface_flattened mesh_surface_flattened
  #                 KIND unit
  #                 SOURCE test/Main.cc test/mesh_surface_flattened.cc
  #                 LINK_LIBS ${mesh_frameworks} mesh_extracted geometry ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

  # note: this requires volumetric extraction, which is only supported by MSTK for now...
  if (ENABLE_MESH_MSTK)
    # Test: column_mesh
    add_amanzi_test(mesh_extracted_column mesh_extracted_column
                  KIND unit
                  SOURCE test/Main.cc test/mesh_extracted_column.cc 
		  LINK_LIBS mesh_extracted mesh_mstk geometry ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # # Test: extracted surface from column_mesh
    # add_amanzi_test(mesh_extracted_surface_cell mesh_extracted_surface_cell
    #               KIND unit
    #               SOURCE test/Main.cc test/mesh_extracted_surface_cell.cc 
    #               LINK_LIBS mesh_extracted mesh_mstk geometry ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})
  endif()          
endif()
  
