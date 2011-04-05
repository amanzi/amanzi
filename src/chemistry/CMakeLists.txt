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

  add_executable(chemistry_tests
		 test/Main.cc
		 test/test_string_tokenizer.cpp
		 test/test_chemistry_exceptions.cpp
		 test/test_species.cpp
		 test/test_ion_exchange_site.cpp
		 test/test_surface_site.cpp
		 test/test_activity_coefficients.cpp
		 test/test_beaker.cpp)
  target_link_libraries(chemistry_tests chemistry ${UnitTest_LIBRARIES})
  add_test(chemistry chemistry_tests)

  if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
    execute_process(COMMAND ${CMAKE_COMMAND} -E
      copy_directory ${CHEM_SOURCE_DIR}/test ${CHEM_BINARY_DIR}/test)
  endif()

endif()

