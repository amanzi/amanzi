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

set(mesh_functions_src_files BoundaryFunction.cc
  UniqueMeshFunction.cc
  DynamicBoundaryFunction.cc
  MaterialMeshFunction.cc
  #MeshPartition.cc
  CompositeVectorFunction.cc
  CompositeVectorFunctionFactory.cc
  )

set(mesh_functions_inc_files MeshFunction.hh
  BoundaryFunction.hh
  UniqueMeshFunction.hh
  DynamicBoundaryFunction.hh
  MaterialMeshFunction.hh
  #MeshPartition.hh
  CompositeVectorFunction.hh
  CompositeVectorFunctionFactory.hh
  )

set(mesh_functions_link_libs mesh data_structures functions error_handling ${Tpetra_LIBRARIES} ${Teuchos_LIBRARIES})

#
# Library: mesh_functions
#
add_amanzi_library(mesh_functions
  SOURCE ${mesh_functions_src_files}
  HEADERS ${mesh_functions_inc_files}
  LINK_LIBS ${mesh_functions_link_libs}
  )

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

    set(mesh_functions_test_link_libs mesh_functions mesh_factory ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES} ${Boost_LIBRARIES})

# these are currently broken
  add_amanzi_test(mesh_functions_boundary_function mesh_functions_boundary_function
                  KIND unit
                  SOURCE test/mesh_functions_boundary_function.cc
                  LINK_LIBS ${mesh_functions_test_link_libs})

  add_amanzi_test(mesh_functions_material_function mesh_functions_material_function
                  KIND unit
                  SOURCE test/mesh_functions_material_function.cc
                  LINK_LIBS ${mesh_functions_test_link_libs})

  add_amanzi_test(mesh_functions_composite_vector_function mesh_functions_composite_vector_function
                  KIND unit
                  SOURCE test/mesh_functions_composite_vector_function.cc
                  LINK_LIBS ${mesh_functions_test_link_libs})
  
endif()
