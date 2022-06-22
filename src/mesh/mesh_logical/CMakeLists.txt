# -*- mode: cmake -*-

#
#  Amanzi
#    Mesh 
#
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)

#
# Define a project name
# After this command the following varaibles are defined
#   MESH_SOURCE_DIR
#   MESH_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MESH_LOGICAL)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})

# External (TPL) directories
include_directories(${Tpetra_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

#
# Library: mesh
#
set(mesh_logical_src MeshLogical.cc
		     MeshLogicalFactory.cc
                     MeshLogicalAudit.cc
		     MeshEmbeddedLogical.cc
                     )

add_amanzi_library(mesh_logical SOURCE ${mesh_logical_src}
				LINK_LIBS mesh geometry ${Tpetra_LIBRARIES})

#
# Install Header files
#
file(GLOB mesh_logical_inc_files "*.hh")
add_install_include_file(${mesh_logical_inc_files})


if (BUILD_TESTS)
   add_definitions("-DMSTK_HAVE_MPI")
   if (WITH_MSTK_2_20rc1_OR_NEWER)
     add_definitions("-DMSTK_2_20rc1_OR_NEWER")
   endif ()

   if (WITH_MSTK_2_21rc1_OR_NEWER)
     add_definitions("-DMSTK_2_21rc1_OR_NEWER")
   endif ()

    # Copy test directory files if an out of source build
    if (NOT (${MESH_LOGICAL_SOURCE_DIR} EQUAL ${MESH_LOGICAL_BINARY_DIR}) )
        file(GLOB DataFiles "${MESH_LOGICAL_SOURCE_DIR}/test/*.xml"
                            "${MESH_LOGICAL_SOURCE_DIR}/test/*.exo")
        file(COPY ${DataFiles} DESTINATION ${MESH_LOGICAL_BINARY_DIR}/test/)
    endif()

   # Include dirs for test
   include_directories(${MESH_LOGICAL_SOURCE_DIR})
   include_directories(${MESH_LOGICAL_SOURCE_DIR}/test)
   include_directories(${MSTK_SOURCE_DIR})
   include_directories(${GEOMETRY_SOURCE_DIR})
   include_directories(${MSTK_INCLUDE_DIRS})

   # Add UnitTest include directories
   include_directories(${UnitTest_INCLUDE_DIRS})

   # Test link libraries
   set(test_link_libs mesh_logical mesh mesh_mstk geometry ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

   # Test: init
   set(test_src_files test/Main.cc test/demo_mesh.cc test/test_mesh_logical.cc)

   add_amanzi_test(mesh_logical test_mesh_logical
                KIND unit
		SOURCE ${test_src_files}
		LINK_LIBS ${test_link_libs})

endif()
  
