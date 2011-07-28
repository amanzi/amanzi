# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(TestManager)
include(InstallManager)

get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)
get_property(GEOCHEM_LIBS GLOBAL PROPERTY GEOCHEM_LIBS)

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
	       test_pitzer.cc
	       test_secondary_species.cc
	       test_aqueous_equilibrium_complex.cc
	       test_mineral.cc
	       test_mineral_kinetics_factory.cc
	       test_kinetic_rate.cc
	       test_beaker.cc)	       
target_link_libraries(chemistry_tests 
                      error_handling ${GEOCHEM_LIBS} ${UnitTest_LIBRARIES})
add_amanzi_test(chemistry chemistry_tests KIND unit)
#
# Install targets
#
file(GLOB chem_inc "${CHEM_INCLUDES_DIR}/*")
add_install_include_file(${chem_inc})

if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
  execute_process(COMMAND ${CMAKE_COMMAND} -E
                  copy_directory ${CHEM_SOURCE_DIR}/tests_unit ${CHEM_BINARY_DIR}/tests_unit)
endif()


