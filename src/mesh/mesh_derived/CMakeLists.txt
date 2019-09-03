#
#  Amanzi
#    Derived Mesh framework
#


# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)

#
# Define a project name
# After this command the following varaibles are defined
#   MESH_DERIVED_SOURCE_DIR
#   MESH_DERIVED_BINARY_DIR
project(MESH_DERIVED)

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_DERIVED_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MSTK_SOURCE_DIR})

# External (TPL) directories
include_directories(${Epetra_INCLUDE_DIRS})

#
# Library: mesh
#
set(mesh_derived_inc_files
    BlockMapUtils.hh
    MeshDerived.hh
)

add_amanzi_library(derived_mesh
                   SOURCE BlockMapUtils.cc
                          MeshDerived.cc
                   HEADERS ${mesh_derived_inc_files}
                   LINK_LIBS mesh geometry ${Epetra_LIBRARIES})

if (BUILD_TESTS)

   include_directories(${UnitTest_INCLUDE_DIRS})

   # Need to copy files for the tests 
   if (NOT (${MSTK_SOURCE_DIR} EQUAL ${MSTK_BINARY_DIR}))
       file(GLOB DataFiles "${MESH_DERIVED_SOURCE_DIR}/test/*.xml")
       file(COPY ${DataFiles} DESTINATION ${MESH_DERIVED_BINARY_DIR}/test/)
   endif()

   add_amanzi_test(mesh_derived_fracture mesh_derived_fracture
                   KIND unit
                   SOURCE test/Main.cc test/mesh_derived_fracture.cc
                   LINK_LIBS mesh_factory mstk_mesh derived_mesh geometry ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

endif()
  
