# -*- mode: cmake -*-


#
#  classes that help with debugging (DBG)
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(LibraryManager)

#
# Define a project name
# After this command the following varaibles are defined
#   DBG_SOURCE_DIR
#   DBG_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(DBG)

include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${Teuchos_INCLUDE_DIRS})


add_amanzi_library(dbg
                   SOURCE DataDebug.cc
                   HEADERS DataDebug.hh
		   LINK_LIBS mesh geometry error_handling ${Teuchos_LIBRARIES}  
                   )
if (BUILD_TESTS)
  
    # Add UnitTest include directoy
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Copy test directory files if an out of source build
    if (NOT (${DBG_SOURCE_DIR} EQUAL ${DBG_BINARY_DIR}))
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${DBG_SOURCE_DIR}/test ${DBG_BINARY_DIR}/test) 
    endif()

    # Add the utils directory to the include paths
    include_directories(${DBG_SOURCE_DIR})

    ADD_AMANZI_TEST(utils test_dbg
                    KIND unit
                    SOURCE test/Main.cc test/test_dbg.cc
                    LINK_LIBS dbg ${UnitTest_LIBRARIES} )

endif()











