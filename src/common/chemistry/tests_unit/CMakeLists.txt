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

add_amanzi_test(chemistry_native chemistry_native
                KIND unit
                SOURCE Main.cc   
                       chemistry_string_tokenizer.cc
                       chemistry_utilities.cc
                       chemistry_exceptions.cc
                       chemistry_species.cc
                       chemistry_surface_site.cc           
                       chemistry_activity_model_factory.cc
                       chemistry_activity_coefficients.cc
                       chemistry_pitzer.cc
                       chemistry_secondary_species.cc
                       chemistry_aqueous_equilibrium_complex.cc
                       chemistry_mineral.cc
                       chemistry_mineral_kinetics_factory.cc
                       chemistry_kinetic_rate.cc
                       chemistry_beaker.cc           
                 LINK_LIBS
                      error_handling atk ${GEOCHEM_LIBS}
                      ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${UnitTest_LIBRARIES})
#
# Install targets
#
file(GLOB chem_inc "${CHEM_INCLUDES_DIR}/*")
add_install_include_file(${chem_inc})

if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
  file(GLOB DataFiles "${CHEM_SOURCE_DIR}/tests_unit/*.bgd"
                      "${CHEM_SOURCE_DIR}/tests_unit/*.xml"
                      "${CHEM_SOURCE_DIR}/tests_unit/*.dat")
  file(COPY ${DataFiles} DESTINATION ${CHEM_BINARY_DIR}/tests_unit/)
endif()


