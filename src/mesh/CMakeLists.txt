# -*- mode: cmake -*-
#
#  Amanzi Mesh
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
# Set files
#
set(mesh_inc_files
  MeshDefs.hh
  CellTopology.hh
  GenerationSpec.hh
  MeshFramework.hh
  Mesh_Helpers_decl.hh
  Mesh_Helpers_impl.hh
  Mesh_Helpers.hh
  MeshAudit_decl.hh
  MeshAudit_impl.hh
  MeshFrameworkAudit.hh
  MeshAudit.hh
  MeshMaps_decl.hh
  MeshMaps_impl.hh
  MeshMaps.hh
  MeshSets.hh
  MeshColumns.hh
  MeshCache.hh
  Mesh.hh
  Mesh_Algorithms.hh
  DomainSet.hh
)

set(mesh_src_files
  GenerationSpec.cc
  MeshFramework.cc
  MeshSets.cc
  MeshColumns.cc
  MeshCache.cc
  Mesh_Algorithms.cc
  DomainSet.cc
)

#
# Library: mesh
#
add_amanzi_library(mesh SOURCE ${mesh_src_files}
                        HEADERS ${mesh_inc_files}
                        LINK_LIBS geometry atk ${Epetra_LIBRARIES})

# Mesh Frameworks
add_subdirectory(mesh_simple)

# MOAB Mesh framework
if (ENABLE_MESH_MOAB)
  add_subdirectory(mesh_moab)
endif(ENABLE_MESH_MOAB)    

# MSTK Mesh framework
if (ENABLE_MESH_MSTK)
  add_subdirectory(mesh_mstk)
endif(ENABLE_MESH_MSTK)    

#add_subdirectory(mesh_extracted)
#add_subdirectory(mesh_logical)

# Factory class
add_subdirectory(mesh_factory)

if (BUILD_TESTS)
    # include mesh frameworks
    include_directories(${MESH_SOURCE_DIR})
    include_directories(${MESH_SIMPLE_SOURCE_DIR})
    if (ENABLE_MESH_MOAB)
      include_directories(${MOAB_INCLUDE_DIR})
      include_directories(${MESH_MOAB_SOURCE_DIR})
      list(APPEND mesh_factory_frameworks mesh_moab)
      add_definitions("-DHAVE_MESH_MOAB")
    endif()
    
    if (ENABLE_MESH_MSTK)
      include_directories(${MSTK_INCLUDE_DIR})
      include_directories(${MESH_MSTK_SOURCE_DIR})
      list(APPEND mesh_factory_frameworks mesh_mstk)
      add_definitions("-DHAVE_MESH_MSTK -DMSTK_HAVE_MPI")
    endif()

    include_directories(${MESH_EXTRACTED_SOURCE_DIR})
    include_directories(${MESH_FACTORY_SOURCE_DIR})
    include_directories(${GEOMETRY_SOURCE_DIR})
    include_directories(${Epetra_INCLUDE_DIRS})
    include_directories(${UnitTest_INCLUDE_DIRS})
    
    # Need to copy files for the tests 
    if (NOT (${MESH_SOURCE_DIR} EQUAL ${MESH_BINARY_DIR}))
      file(GLOB DataFiles "${MESH_SOURCE_DIR}/test/*.exo"
        "${MESH_SOURCE_DIR}/test/*.par.*"
        "${MESH_SOURCE_DIR}/test/*.xml"
        "${MESH_SOURCE_DIR}/test/*.txt" )

      file(COPY ${DataFiles} DESTINATION ${MESH_BINARY_DIR}/test/)
    endif()

    # Test link libraries
    set(mesh_test_link_libs geometry mesh mesh_simple mesh_factory ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

    # MOAB Mesh framework
    if (ENABLE_MESH_MOAB)
      list(APPEND mesh_test_link_libs mesh_moab)
    endif()    
    
    # MSTK Mesh framework
    if (ENABLE_MESH_MSTK)
      list(APPEND mesh_test_link_libs mesh_mstk)
    endif()    

    # Test: tests framework geometry on a variety of meshes
    add_amanzi_test(mesh_framework_geometry mesh_framework_geometry
                    KIND unit
                    SOURCE test/Main.cc test/mesh_framework_geometry.cc
                    LINK_LIBS ${mesh_test_link_libs})
    add_amanzi_test(mesh_framework_geometry_np2 mesh_framework_geometry NPROCS 2 KIND unit)
    add_amanzi_test(mesh_framework_geometry_np4 mesh_framework_geometry NPROCS 4 KIND unit)

    # Test: tests MeshCache geometry without a cache on a variety of meshes
    add_amanzi_test(mesh_geometry mesh_geometry
                    KIND unit
                    SOURCE test/Main.cc test/mesh_geometry.cc
                    LINK_LIBS ${mesh_test_link_libs})
    add_amanzi_test(mesh_geometry_np2 mesh_geometry NPROCS 2 KIND unit)
    add_amanzi_test(mesh_geometry_np4 mesh_geometry NPROCS 4 KIND unit)

    # Test: tests MeshCace geometry with a cache and no framework on a variety of meshes
    add_amanzi_test(mesh_cache_geometry mesh_cache_geometry
                    KIND unit
                    SOURCE test/Main.cc test/mesh_cache_geometry.cc
                    LINK_LIBS ${mesh_test_link_libs})
    add_amanzi_test(mesh_cache_geometry_np2 mesh_cache_geometry NPROCS 2 KIND unit)
    add_amanzi_test(mesh_cache_geometry_np4 mesh_cache_geometry NPROCS 4 KIND unit)
    
    # Test: mesh sets
    #   Tests the use of mesh sets on the above meshes and frameworks.
    add_amanzi_test(mesh_sets mesh_sets
                    KIND unit
                    SOURCE test/Main.cc test/mesh_sets.cc
                    LINK_LIBS ${mesh_test_link_libs})
    add_amanzi_test(mesh_sets_np2 mesh_sets NPROCS 2 KIND unit)
    add_amanzi_test(mesh_sets_np4 mesh_sets NPROCS 4 KIND unit)
    
    # Test: mesh deformation
    #   A basic test of deforming volumetric and surface meshes and
    #   making sure geometric quantities get updated correctly.
    add_amanzi_test(mesh_deform mesh_deform
                    KIND unit
        	    SOURCE test/Main.cc test/mesh_deform.cc
        	    LINK_LIBS ${mesh_test_link_libs})

    # Test: mesh columns
    #   Tests the buidl_columns() capability for semi-structured meshes.              
    add_amanzi_test(mesh_build_columns mesh_build_columns
                    KIND unit
        	    SOURCE test/Main.cc test/mesh_build_columns.cc
        	    LINK_LIBS ${mesh_test_link_libs})
    add_amanzi_test(mesh_build_columns_np2 mesh_build_columns NPROCS 2 KIND unit)

    # Test: surface extracted mesh
    #   Tests a surface mesh generated using the mesh framework
    #   extraction constructor on FACEs with the flatten parameter,
    #   currently only supported by MSTK.  This extracts a surface
    #   mesh and checks MeshAudit, geometry, sets, and maps.
    #
    # NOTE: Currently, the MeshExtractedManifold tests in here are commented out!
    add_amanzi_test(mesh_surface_extraction mesh_surface_extraction
                    KIND unit
                    SOURCE test/Main.cc test/mesh_surface_extraction.cc
                    LINK_LIBS ${mesh_test_link_libs})
    add_amanzi_test(mesh_surface_extraction_np2 mesh_surface_extraction NPROCS 2 KIND unit)
    add_amanzi_test(mesh_surface_extraction_np4 mesh_surface_extraction NPROCS 4 KIND unit)
    
    # Test: volume extracted mesh
    #   Tests a volumetric mesh generated using the mesh framework
    #   extraction constructor on CELLs, currently only supported by
    #   MSTK.  This extracts a set of cells, forms a mesh, and
    #   checks MeshAudit, geometry, sets, and maps.
    add_amanzi_test(mesh_volume_extraction mesh_volume_extraction
                    KIND unit
                    SOURCE test/Main.cc test/mesh_volume_extraction.cc
                    LINK_LIBS ${mesh_test_link_libs})
    
    # Test: fracture extracted mesh
    #   Tests a fracture mesh generated using the mesh framework
    #   extraction constructor on FACES but not flattened, currently
    #   only supported by MSTK.  This extracts a set of faces and
    #   should recognize the shared edges between two intersecting
    #   planes.
    # add_amanzi_test(mesh_fracture_extraction mesh_fracture_extraction
    #                 KIND unit
    #                 SOURCE test/Main.cc test/mesh_fracture_extraction.cc
    #                 LINK_LIBS ${mesh_test_link_libs})

    # Test: mesh KD search
    add_amanzi_test(mesh_nanoflann mesh_nanoflann
                    KIND unit
        	    SOURCE test/Main.cc test/mesh_nanoflann.cc
        	    LINK_LIBS ${mesh_test_link_libs})
endif()
  
