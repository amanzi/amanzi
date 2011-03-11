# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(AddParallelTest)

#
# Define a project name
# After this command the following varaibles are defined
#   CHEM_SOURCE_DIR
#   CHEM_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(CHEM)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_BASE_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})

# Chemistry includes directory
# The set will define the variable in the chemistry scope
# The global property will allow other directories (MPC) to
# access this directory name through a variable.
set(CHEM_INCLUDES_DIR "${CHEM_SOURCE_DIR}/includes")
set_property(GLOBAL PROPERTY CHEM_INCLUDES_DIR "${CHEM_INCLUDES_DIR}")
include_directories(${CHEM_INCLUDES_DIR})



add_subdirectory(utility)
add_subdirectory(solvers)
add_subdirectory(base_chemistry)
add_subdirectory(reactions)
add_subdirectory(process_kernel)
add_subdirectory(test_drivers)


if (BUILD_TESTS)
  # Add the UnitTest includes
  include_directories(${UnitTest_INCLUDE_DIRS})

  # External (TPL) includes
  include_directories(${Teuchos_INCLUDE_DIRS})

  add_executable(test_activity_coefficients test/Main.cc test/test_activity_coefficients.cpp)
  add_test(chemistry_activity_coefficients test_activity_coefficients)
  target_link_libraries(test_activity_coefficients chemistry ${UnitTest_LIBRARIES})

  add_executable(test_beaker test/Main.cc test/test_beaker.cpp)
  add_test(chemistry_beaker test_beaker)
  target_link_libraries(test_beaker chemistry ${UnitTest_LIBRARIES})


    # Test: chemistry_unit
    # add_executable(test_chemistry_unit test/Main.cc test/test_example.cpp)
    # target_link_libraries(test_chemistry_unit chemistry ${UnitTest_LIBRARIES} ${Teuchos_LIBRARIES})
    # add_test(chemistry_unit test_chemistry_unit)

    # Test: chemistry
    #add_executable(test_chemistry test/Main.cc test/test_driver.cpp)
    #add_test(chemistry test_chemistry)
    #target_link_libraries(test_chemistry state simple_mesh mpc chemistry transport flow ${UnitTest_LIBRARIES})


    if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
            copy_directory ${CHEM_SOURCE_DIR}/test ${CHEM_BINARY_DIR}/test) 
    endif()

endif()

