# -*- mode: cmake -*-
#
#  AmanziGeometry
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

#
# Define a project name
# After this command the following varaibles are defined
#   GEOMETRY_SOURCE_DIR
#   GEOMETRY_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(GEOMETRY)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})  

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})

#
# geometry library
#
file(GLOB geometry_inc_files "*.hh")
file(GLOB geometry_src_files "*.cc")

add_amanzi_library(geometry 
  SOURCE ${geometry_src_files}
  HEADERS ${geometry_inc_files}
  LINK_LIBS functions error_handling atk ${Teuchos_LIBRARIES} )


#
# Unit tests
#
if (BUILD_TESTS)
   include_directories(${GEOMETRY_SOURCE_DIR})
   include_directories(${UnitTest_INCLUDE_DIRS})

   # Copy test directory files if an out of source build
   if (NOT (${GEOMETRY_SOURCE_DIR} EQUAL ${GEOMETRY_BINARY_DIR}) )
       file(GLOB DataFiles "${GEOMETRY_SOURCE_DIR}/test/*.xml"
                           "${GEOMETRY_SOURCE_DIR}/test/*.txt"
                           "${GEOMETRY_SOURCE_DIR}/test/*.exo")
       file(COPY ${DataFiles} DESTINATION ${GEOMETRY_BINARY_DIR}/test/)
   endif()

   # Test: test point creation and deletion
   add_amanzi_test(geometry_points geometry_points
                   KIND unit
		   SOURCE test/Main.cc test/test_points.cc
		   LINK_LIBS geometry ${UnitTest_LIBRARIES})

   # Test: test geometric operators
   add_amanzi_test(geometry_operators geometry_operations
                   KIND unit
                   SOURCE test/Main.cc 
                          test/test_geometric_ops.cc
                   LINK_LIBS geometry ${UnitTest_LIBRARIES})

   # Test: test region creation
   add_amanzi_test(geometry_region geometry_region
                   KIND unit
                   SOURCE test/Main.cc 
                          test/test_box_region.cc
                          test/test_plane_region.cc
                          test/test_labeledset.cc
                          test/test_enumeratedset.cc
                          test/test_colorfunc_region.cc
			  test/test_polygon_region.cc
			  test/test_cylinder_region.cc
                   LINK_LIBS geometry ${UnitTest_LIBRARIES})

   # Test: test geometric model creation
   add_amanzi_test(geometry_geometric_model geometry_geometric_model
                   KIND unit
                   SOURCE test/Main.cc
                          test/test_geometric_model.cc
                   LINK_LIBS geometry
                          ${UnitTest_LIBRARIES} ${Teuchos_LIBRARIES})  

endif()

