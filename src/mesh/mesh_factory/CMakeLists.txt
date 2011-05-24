# -*- mode: cmake -*-------------------------------------------
# file: CMakeLists.txt
#
# Amanzi
#   Mesh Factory
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created March 11, 2011 by William A. Perkins
# Last Change: Fri Apr 22 09:49:51 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
# -------------------------------------------------------------
# Amanzi module, include files found in AMANZI_MODULE_PATH
include(TestManager)
include(PrintVariable)


#
# Define a project name
# After this command the following varaibles are defined
#   MESH_FACTORY_SOURCE_DIR
#   MESH_FACTORY_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MESH_FACTORY)

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${STK_INCLUDE_DIRS})
include_directories(${ExodusII_INCLUDE_DIRS})

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_SIMPLE_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})

#
# Library: mesh_factory
#

set(mesh_factory_files MeshFileType.cc MeshFramework.cc MeshFactory.cc FrameworkTraits.cc) 
add_library(mesh_factory ${mesh_factory_files})
target_link_libraries(mesh_factory error_handling simple_mesh)

if (ENABLE_STK_Mesh)
    include_directories(${STK_INCLUDE_DIR} ${STK_SOURCE_DIR})
    target_link_libraries(mesh_factory stk_mesh exodus mesh_data)
    add_definitions("-DHAVE_STK_MESH")
endif()

if (ENABLE_MOAB_Mesh)
    include_directories(${MOAB_INCLUDE_DIR} ${MOAB_SOURCE_DIR})
    target_link_libraries(mesh_factory moab_mesh)
    add_definitions("-DHAVE_MOAB_MESH")
endif()

if (ENABLE_MSTK_Mesh)
    include_directories(${MSTK_INCLUDE_DIR} ${MSTK_SOURCE_DIR})
    target_link_libraries(mesh_factory mstk_mesh)
    add_definitions("-DHAVE_MSTK_MESH -DMSTK_HAVE_MPI")
endif()

if (BUILD_TESTS)

    # some test file locations defined as macros

    # files in the exodus test directory always seem to be available

    set(EXODUS_FILE_DIR "${EXO_BINARY_DIR}/test_files")
    set(test_file_flags "-DEXODUS_TEST_FILE=\"\\\"${EXODUS_FILE_DIR}/hex_4x4x4_ss.exo\\\"\" ")
    set(test_file_flags "${test_file_flags} -DNEMESIS_TEST_FILE=\"\\\"${EXODUS_FILE_DIR}/split1/hex_11x11x11_ss.par\\\"\" ")
    set(test_file_flags "${test_file_flags} -DBOGUS_TEST_FILE=\"\\\"${EXODUS_FILE_DIR}/bogus.exo\\\"\" ")

    # If the MOAB framework is not available, MOAB_SOURCE_DIR is not
    # defined.  We still need a MOAB file to check, though.  Wouldn't
    # it be better to have a single place where we put all
    # example/test meshes?
    if(ENABLE_MOAB_Mesh)
        set(test_file_flags "${test_file_flags} -DMOAB_TEST_FILE=\"\\\"${MOAB_SOURCE_DIR}/test/hex_4x4x4_ss_4P.h5m\\\"\" ")
    else()
        set(test_file_flags "${test_file_flags} -DMOAB_TEST_FILE=\"\\\"${AMANZI_SOURCE_DIR}/src/moab_mesh/test/hex_4x4x4_ss_4P.h5m\\\"\" ")
    endif()
        
    set_source_files_properties(test/test_mesh_file.cc PROPERTIES 
                                COMPILE_FLAGS ${test_file_flags})

    set_source_files_properties(test/test_mesh_factory.cc PROPERTIES 
                                COMPILE_FLAGS ${test_file_flags})
   
   # Add UnitTest include directories
   include_directories(${UnitTest_INCLUDE_DIRS})

   # Test: test_mesh_factory
   add_executable(test_mesh_factory 
     test/Main.cc
     test/test_mesh_file.cc 
     test/test_mesh_framework.cc
     test/test_mesh_factory.cc)
   target_link_libraries(test_mesh_factory mesh_factory ${UnitTest_LIBRARIES} ${Boost_LIBRARIES})
   add_amanzi_test(mesh_factory test_mesh_factory KIND int)
   add_amanzi_test(mesh_factory_parallel test_mesh_factory KIND int)

endif()
   
