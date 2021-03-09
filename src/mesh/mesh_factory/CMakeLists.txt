# -*- mode: cmake -*-------------------------------------------
# file: CMakeLists.txt
#
# Amanzi
#   Mesh Factory
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created March 11, 2011 by William A. Perkins
# Last Change: Wed Sep 28 09:07:13 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
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
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${STK_INCLUDE_DIRS})
include_directories(${SEACAS_INCLUDE_DIRS})

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_EXTRACTED_SOURCE_DIR})
include_directories(${MESH_SIMPLE_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})


#
# Library: mesh_factory
#
set(mesh_factory_files 
    FileFormat.cc
    MeshFrameworkTraits.cc
    MeshFrameworkFactory.cc) 

set(mesh_factory_inc_files
    MeshException.hh
    FileFormat.hh
    MeshFrameworkTraits.hh
    MeshFrameworkFactory.hh)


list(APPEND mesh_factory_frameworks simple_mesh)

if (ENABLE_STK_Mesh)
    include_directories(${STK_INCLUDE_DIR} ${STK_SOURCE_DIR})
    list(APPEND mesh_factory_frameworks stk_mesh exodus_reader mesh_data)
    add_definitions("-DHAVE_STK_MESH")
endif()

if (ENABLE_MOAB_Mesh)
    include_directories(${MOAB_INCLUDE_DIR} ${MOAB_SOURCE_DIR})
    list(APPEND mesh_factory_frameworks moab_mesh)
    add_definitions("-DHAVE_MOAB_MESH")
endif()

if (ENABLE_MSTK_Mesh)
    include_directories(${MSTK_INCLUDE_DIR} ${MSTK_SOURCE_DIR})
    list(APPEND mesh_factory_frameworks mstk_mesh)
    add_definitions("-DHAVE_MSTK_MESH -DMSTK_HAVE_MPI")
endif()

add_amanzi_library(mesh_factory
                   SOURCE ${mesh_factory_files} HEADERS ${mesh_factory_inc_files}
                   LINK_LIBS ${mesh_factory_frameworks} ${Boost_LIBRARIES})
#                   LINK_LIBS ${mesh_factory_frameworks} mesh_extracted ${Boost_LIBRARIES})

if (BUILD_TESTS)

    # Need to copy files for the tests 
    if (NOT (${MESH_FACTORY_SOURCE_DIR} EQUAL ${MESH_FACTORY_BINARY_DIR}))
        file(GLOB DataFiles "${MESH_FACTORY_SOURCE_DIR}/test/*.xml"
                            "${MESH_FACTORY_SOURCE_DIR}/test/*.txt"
                            "${MESH_FACTORY_SOURCE_DIR}/test/*.h5m"
                            "${MESH_FACTORY_SOURCE_DIR}/test/*.exo"
                            "${MESH_FACTORY_SOURCE_DIR}/test/*.par.*")
        file(COPY ${DataFiles} DESTINATION ${MESH_FACTORY_BINARY_DIR}/test/)
    endif()

   # Add UnitTest include directories
   include_directories(${UnitTest_INCLUDE_DIRS})

   # Test: test_mesh_factory
   add_amanzi_test(mesh_factory test_mesh_factory
                   KIND unit
		   SOURCE
                         test/Main.cc
                         test/test_mesh_file.cc 
                         test/test_mesh_factory.cc
                   LINK_LIBS
		         mesh_factory
			 ${UnitTest_LIBRARIES}
			 ${Boost_LIBRARIES}
                         ${Epetra_LIBRARIES}
			 ${HDF5_LIBRARIES})

   # Do not call add_amanzi_test again with SOURCE and LINK_LIBS since the
   # previous call carries instructions to build the executable
   add_amanzi_test(mesh_factory_parallel test_mesh_factory NPROCS 2 KIND uint)

endif()
   
