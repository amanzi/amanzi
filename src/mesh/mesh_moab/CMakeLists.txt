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
project(MOAB)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})

# Need this define. Errors from MOAB include files 
# about MPI_COOM_WORLD. --lpritch
add_definitions("-DUSE_MPI")

#
# Library: moab_mesh
#
# The MOAB_LIBRARIES variable does not have HDF5
# libraries listed. Ticket #335
file(GLOB moab_mesh_files "*.cc") 
file(GLOB moab_inc_files "*.hh")
add_amanzi_library(moab_mesh
                   SOURCE ${moab_mesh_files} HEADERS ${moab_inc_files}
                   LINK_LIBS ${MOAB_LIBRARIES} ${Epetra_LIBRARIES} ${Teuchos_LIBRARIES}) 

if (BUILD_TESTS)
    
    # Add UnitTest include
#    include_directories(${UnitTest_INCLUDE_DIRS})

    # Need to copy files for the tests 
#    if (NOT (${MOAB_SOURCE_DIR} EQUAL ${MOAB_BINARY_DIR}))
#        execute_process(COMMAND ${CMAKE_COMMAND} -E 
#          copy_directory ${MOAB_SOURCE_DIR}/test ${MOAB_BINARY_DIR}/test) 
#    endif()


    # Test: moab_Mesh
#    add_executable(test_moab_mesh 
#                   test/Main.cc 
#                   test/test_hex_2x2x2.cc
#                   test/test_hex_3x3x2.cc
#                   test/test_hex_4x4x4.cc)
#    target_link_libraries(test_moab_mesh moab_mesh ${UnitTest_LIBRARIES})
#    add_amanzi_test(moab_mesh test_moab_mesh KIND unit)
#
    # Test: moab_mesh_parallel
#    add_executable(test_parallel_moab_mesh
#        test/Main.cc
#        test/test_hex_4x4x4_4P.cc)
#    target_link_libraries(test_parallel_moab_mesh moab_mesh ${UnitTest_LIBRARIES})
#    add_amanzi_test(moab_mesh_parallel test_parallel_moab_mesh NPROCS 4 KIND int)

endif()

