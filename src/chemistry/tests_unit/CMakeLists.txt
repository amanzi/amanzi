# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(TestManager)

# Chemistry includes
include_directories(${CHEM_INCLUDES_DIR})

# Add the UnitTest includes
include_directories(${UnitTest_INCLUDE_DIRS})

add_executable(chemistry_tests
               main.cc
	       test_string_tokenizer.cc
	       test_chemistry_exceptions.cc
	       test_species.cc
	       test_ion_exchange_site.cc
	       test_surface_site.cc
	       test_activity_coefficients.cc
	       test_secondary_species.cc
	       test_aqueous_equilibrium_complex.cc
	       test_mineral.cc
	       test_mineral_kinetics_factory.cc
	       test_kinetic_rate.cc
	       test_beaker.cc)

target_link_libraries(chemistry_tests chemistry 
                      error_handling ${UnitTest_LIBRARIES})
add_amanzi_test(chemistry chemistry_tests KIND unit)


# Trilinos related (TPL) includes
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})

add_executable(chemistry_pk_tests
               main.cc
	       test_chemistry_pk.cc)

target_link_libraries(chemistry_pk_tests chemistry state 
                      simple_mesh error_handling exodus 
                      ${UnitTest_LIBRARIES} 
		      ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES})
add_amanzi_test(chemistry_pk chemistry_pk_tests KIND unit)

if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
  execute_process(COMMAND ${CMAKE_COMMAND} -E
                  copy_directory ${CHEM_SOURCE_DIR}/tests_unit ${CHEM_BINARY_DIR}/tests_unit)
endif()


