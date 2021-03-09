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
project(MESH)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})

# External (TPL) directories
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

#
# Library: mesh
#
add_amanzi_library(mesh SOURCE GenerationSpec.cc
                               MeshFramework.cc
                               #DomainSet.cc
                        LINK_LIBS geometry atk ${Epetra_LIBRARIES})

#
# Library: mesh_audit
#
add_amanzi_library(mesh_audit SOURCE MeshAudit.cc 
	                      LINK_LIBS mesh error_handling ${Epetra_LIBRARIES})

#
# Install Header files
#
file(GLOB mesh_inc_files "*.hh")
add_install_include_file(${mesh_inc_files})


# Mesh Frameworks
add_subdirectory(mesh_simple)

# MOAB Mesh framework
if (ENABLE_MOAB_Mesh)
#  add_subdirectory(mesh_moab)
endif(ENABLE_MOAB_Mesh)    

# MSTK Mesh framework
if (ENABLE_MSTK_Mesh)
  add_subdirectory(mesh_mstk)
endif(ENABLE_MSTK_Mesh)    

#add_subdirectory(mesh_extracted)
#add_subdirectory(mesh_logical)

# Factory class
add_subdirectory(mesh_factory)


if (BUILD_TESTS)
    # Include dirs for test
    include_directories(${MESH_SOURCE_DIR})
    include_directories(${MESH_FACTORY_SOURCE_DIR})
    include_directories(${GEOMETRY_SOURCE_DIR})
    include_directories(${Epetra_INCLUDE_DIRS})
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Need to copy files for the tests 
    if (NOT (${MESH_SOURCE_DIR} EQUAL ${MESH_BINARY_DIR}))
      # file(GLOB DataFiles "${MESH_SOURCE_DIR}/test/*.xml")
      # file(COPY ${DataFiles} DESTINATION ${MESH_BINARY_DIR}/test/)
    endif()

    # Test link libraries
    set(test_link_libs simple_mesh geometry mesh_factory ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # MOAB Mesh framework
    if (ENABLE_MOAB_Mesh)
      list(APPEND test_link_libs moab_mesh)
    endif()    
    
    # MSTK Mesh framework
    if (ENABLE_MSTK_Mesh)
      list(APPEND test_link_libs mstk_mesh)
    endif()    

    # Test: mesh geometry
    add_amanzi_test(mesh_geometry mesh_geometry
                    KIND unit
                    SOURCE test/Main.cc test/test_mesh_geometry.cc
                    LINK_LIBS ${test_link_libs})

    # Test: mesh deformation
    # add_amanzi_test(mesh_deform mesh_deform
    #                 KIND unit
    #     	    SOURCE test/Main.cc test/test_deform.cc
    #     	    LINK_LIBS ${test_link_libs})

    # Test: mesh columns
    # add_amanzi_test(mesh_columns mesh_columns
    #                 KIND unit
    #     	    SOURCE test/Main.cc test/test_columns.cc
    #     	    LINK_LIBS ${test_link_libs})
    # add_amanzi_test(mesh_columns_np2 mesh_columns NPROCS 2 KIND uint)

    # Test: mesh KD search
    # add_amanzi_test(mesh_nanoflann mesh_nanoflann
    #                 KIND unit
    #     	    SOURCE test/Main.cc test/test_nanoflann.cc
    #     	    LINK_LIBS ${test_link_libs})
endif()
  
