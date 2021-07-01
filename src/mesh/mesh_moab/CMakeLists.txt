# -*- mode: cmake -*-

#
#  Amanzi
#    MOAB Mesh frameworks
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   MOAB_SOURCE_DIR
#   MOAB_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MESH_MOAB)

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})

# Need this define. Errors from MOAB include files 
# about MPI_COOM_WORLD. --lpritch
add_definitions("-DUSE_MPI")

#
# Library: mesh_moab
#
# The MOAB_LIBRARIES variable does not have HDF5
# libraries listed. Ticket #335
file(GLOB mesh_moab_src_files "*.cc") 
file(GLOB mesh_moab_inc_files "*.hh")
add_amanzi_library(mesh_moab
                   SOURCE ${mesh_moab_src_files} HEADERS ${mesh_moab_inc_files}
                   LINK_LIBS geometry mesh atk
                             ${MOAB_LIBRARIES} ${Epetra_LIBRARIES} ${Teuchos_LIBRARIES}) 

if (BUILD_TESTS)
    # Add UnitTest include
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Need to copy files for the tests 
    if (NOT (${MESH_MOAB_SOURCE_DIR} EQUAL ${MESH_MOAB_BINARY_DIR}))
      file(GLOB DataFiles "${MESH_MOAB_SOURCE_DIR}/test/*.exo"
                          "${MESH_MOAB_SOURCE_DIR}/test/*.h5m")
      file(COPY ${DataFiles} DESTINATION ${MESH_MOAB_BINARY_DIR}/test/)
    endif()

    # Test: serial
    add_amanzi_test(mesh_moab_serial mesh_moab_serial
      KIND unit
      SOURCE  test/Main.cc
	      test/mesh_moab_hex_1x1x1.cc
              test/mesh_moab_hex_2x2x1.cc
              test/mesh_moab_hex_3x3x3.cc
      LINK_LIBS mesh_moab mesh_audit ${UnitTest_LIBRARIES})
 
    # Test: parallel
    add_amanzi_test(mesh_moab_parallel mesh_moab_parallel
      KIND unit
      NPROCS 4 
      SOURCE  test/Main.cc
              test/mesh_moab_hex_3x3x3_4P.cc
      LINK_LIBS mesh_moab mesh_audit ${UnitTest_LIBRARIES})

endif()

