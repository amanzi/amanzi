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
#   GEOMETRY_SOURCE_DIR
#   GEOMETRY_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(GEOMETRY)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR} ${FUNCS_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})

#
# geometry library
#
file (GLOB geometry_source_files "*.cc")
add_library(geometry ${geometry_source_files})
target_link_libraries(geometry error_handling functions ${Teuchos_LIBRARIES})

#
# Install Targets
#
add_install_library(geometry)
file(GLOB _inc_files "*.h*")
add_install_include_file(${_inc_files})


#
# Unit tests
#

if (BUILD_TESTS)
    
   # Add UnitTest ExodusII include directories
   include_directories(${UnitTest_INCLUDE_DIRS})

   # Need to copy files for the tests 
   if (NOT (${GEOMETRY_SOURCE_DIR} EQUAL ${GEOMETRY_BINARY_DIR}))
       execute_process(COMMAND ${CMAKE_COMMAND} -E 
         copy_directory ${GEOMETRY_SOURCE_DIR}/test ${GEOMETRY_BINARY_DIR}/test) 
   endif()

   # Test: test geometric operators

   add_executable(test_geomops
                  test/Main.cc
		  test/test_geometric_ops.cc)
   target_link_libraries(test_geomops geometry ${Teuchos_LIBRARIES} ${Trilinos_LIBRARIES} ${UnitTest_LIBRARIES})

   # Test: test region creation
   add_executable(test_region 
                  test/Main.cc 
                  test/test_box_region.cc
                  test/test_plane_region.cc
                  test/test_labeledset.cc
		  test/test_colorfunc_region.cc
		  )
   target_link_libraries(test_region geometry ${Teuchos_LIBRARIES} ${Trilinos_LIBRARIES} ${UnitTest_LIBRARIES})

   # Test: test geometric model creation
   add_executable(test_geometric_model
                  test/Main.cc 
                  test/test_geometric_model.cc)
   target_link_libraries(test_geometric_model geometry ${Teuchos_LIBRARIES} ${Trilinos_LIBRARIES} ${UnitTest_LIBRARIES})


   add_amanzi_test(geometry test_region test_geometric_model test_geomops KIND unit)

    # Test: test geometric model creation from regions
#    add_executable(test_geometric_model
#        test/Main.cc
#        test/test_geometric_model.cc)
#    target_link_libraries(test_geometric_model geometry ${UnitTest_LIBRARIES})
#    add_amanzi_test(geometry test_geometric_model KIND int)

endif()

