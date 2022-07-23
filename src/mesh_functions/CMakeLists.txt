# -*- mode: cmake -*-

#
#  Amanzi
#    Mesh functions
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(LibraryManager)

#
# Define a project name
# After this command the following variables are defined
#   XXX_SOURCE_DIR
#   XXX_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MFUNCS)

# Define include directories to build any binary or library
# in this directory. 

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})

#
# Library: mesh_functions
#
add_amanzi_library(mesh_functions
                   SOURCE CompositeVectorFunction.cc
                          CompositeVectorFunctionFactory.cc
                          UniqueMeshFunction.cc
                          BoundaryFunction.cc
                          ColumnMeshFunction.cc
                          DynamicBoundaryFunction.cc
                          MeshPartition.cc
                          MaterialMeshFunction.cc
                   HEADERS BoundaryFunction.hh
                           ColumnMeshFunction.hh
                           DynamicBoundaryFunction.hh
                           CompositeVectorFunction.hh
                           CompositeVectorFunctionFactory.hh
                           MeshPartition.hh
                           MeshFunction.hh
                           UniqueMeshFunction.hh
                           MaterialMeshFunction.hh
                   LINK_LIBS mesh data_structures error_handling functions ${Teuchos_LIBRARIES})

# Define all tests at the of the file. BUILD_TESTS gates test building.
if (BUILD_TESTS)
    
  # Add UnitTest includes
  include_directories(${UnitTest_INCLUDE_DIRS})

  # Add include directory. Need to remove this required path. -- lpritch
  include_directories(${MFUNCS_SOURCE_DIR})
  include_directories(${FUNCS_SOURCE_DIR})
  include_directories(${MESH_FACTORY_SOURCE_DIR})

  # Copy test subdirectory for out of source builds
  if (NOT (MFUNCS_SOURCE_DIR STREQUAL MFUNCS_BINARY_DIR))
      # file(GLOB DataFiles "${MFUNCS_SOURCE_DIR}/test/*.xml")
      # file(COPY ${DataFiles} DESTINATION ${MFUNCS_BINARY_DIR}/test/)
  endif()

# these are currently broken
  add_amanzi_test(boundary_function boundary_function
                  KIND unit
                  SOURCE test/boundary_function.cc
                  LINK_LIBS mesh_functions mesh_factory data_structures
                            ${UnitTest_LIBRARIES} ${NOX_LIBRARIES} ${HDF5_LIBRARIES})

  add_amanzi_test(material_function material_function
                  KIND unit
                  SOURCE test/material_function.cc
                  LINK_LIBS mesh_functions mesh_factory data_structures
                            ${UnitTest_LIBRARIES} ${NOX_LIBRARIES} ${HDF5_LIBRARIES})

  add_amanzi_test(cv_function cv_function
                  KIND unit
                  SOURCE test/composite_vector_function.cc
                  LINK_LIBS  mesh_functions mesh_factory data_structures 
                             ${UnitTest_LIBRARIES} ${Boost_LIBRARIES} 
                             ${MSTK_LIBRARIES} ${NOX_LIBRARIES} ${HDF5_LIBRARIES})
endif()
