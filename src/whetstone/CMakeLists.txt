# -*- mode: cmake -*-

#
#  Amanzi
#   Discretization library
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)
include(LibraryManager)

#
# Define a project name
# After this command the following varaibles are defined
#   WHETSTONE_SOURCE_DIR
#   WHETSTONE_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(WHETSTONE)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_BASE_SOURCE_DIR})
include_directories(${MESH_DATA_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

# Need to move this to the root level -- lpritch
add_definitions("-DUSE_MPI" "-DMSTK_HAVE_MPI")

#
# Library: whetstone
#
file(GLOB whetstone_inc_files "*.hpp")
add_amanzi_library(whetstone
                   SOURCE mfd3d.cpp mfd3d_vag.cpp tensor.cpp HEADERS ${whetstone_inc_files}
                   LINK_LIBS geometry mesh mesh_data error_handling ${Epetra_LIBRARIES})


if ( BUILD_TESTS )
    # Add UnitTest includes
    include_directories(${UnitTest_INCLUDE_DIRS})
    include_directories(${MESH_SIMPLE_SOURCE_DIR})

    set(amanzi_libs geometry mesh mesh_audit simple_mesh state whetstone)

    # Test 1: inverse mass matrix in darcy flow
    add_amanzi_test(whetstone_darcy whetstone_darcy
                    KIND unit
                    SOURCE test/Main.cc test/whetstone_darcy.cpp
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
endif()


