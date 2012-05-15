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

#
# Define a project name
# After this command the following varaibles are defined
#   ATK_SOURCE_DIR
#   ATK_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(ATK)

include_directories(${Boost_INCLUDE_DIRS})

set(atk_inc_files Utils.hh
                  Timer.hh
                  TimerManager.hh
		  time_step_manager.hh
)

include_directories(${Boost_INCLUDE_DIRS})

add_library(atk Timer.cc
                TimerManager.cc
		time_step_manager.cc
)

add_install_library(atk)

#
# Install Targets
#
add_install_include_file(${atk_inc_files})


if (BUILD_TESTS)
  
    # Add UnitTest include directoy
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Copy test directory files if an out of source build
    if (NOT (${TIME_INTEGRATION_SOURCE_DIR} EQUAL ${TIME_INTEGRATION_BINARY_DIR}))
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${TIME_INTEGRATION_SOURCE_DIR}/test ${TIME_INTEGRATION_BINARY_DIR}/test) 
    endif()

    # Add the utils directory to the include paths
    include_directories(${ATK_SOURCE_DIR})

    add_executable(test_utils test/Main.cc test/test_time_step_manager.cc)    
    target_link_libraries(test_utils atk ${UnitTest_LIBRARIES})
    add_amanzi_test(utils test_utils KIND unit)
endif()











