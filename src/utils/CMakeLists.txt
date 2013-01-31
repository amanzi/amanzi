# -*- mode: cmake -*-


#
#  Amanzi Tool Kit (AKT)
#    General utilities
#    This directory has become a dumping ground for things that do not
#    have a well-defined home. Will go away or get renamed. For now we
#    need to keep it around.
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(LibraryManager)

#
# Define a project name
# After this command the following varaibles are defined
#   ATK_SOURCE_DIR
#   ATK_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(ATK)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})


add_amanzi_library(atk
                   SOURCE Timer.cc TimerManager.cc TimeStepManager.cc XMLParameterListWriter.cc
                   HEADERS Utils.hh Timer.hh TimerManager.hh TimeStepManager.hh XMLParameterListWriter.hh
                   )
if (BUILD_TESTS)
  
    # Add UnitTest include directoy
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Copy test directory files if an out of source build
    if (NOT (${ATK_SOURCE_DIR} EQUAL ${ATK_BINARY_DIR}))
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${ATK_SOURCE_DIR}/test ${ATK_BINARY_DIR}/test) 
    endif()

    # Add the utils directory to the include paths
    include_directories(${ATK_SOURCE_DIR})

    ADD_AMANZI_TEST(utils test_utils
                    KIND unit
                    SOURCE test/Main.cc test/test_time_step_manager.cc
                    LINK_LIBS atk ${UnitTest_LIBRARIES} )

endif()











