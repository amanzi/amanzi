# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)
include(InstallManager)
include(LibraryManager)

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
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})

# Chemistry includes directory
# The set will define the variable in the chemistry scope
# The global property will allow other directories (MPC) to
# access this directory name through a variable.

set(CHEM_INCLUDES_DIR "${CHEM_SOURCE_DIR}/utility"
                      "${CHEM_SOURCE_DIR}/base_chemistry"
                      "${CHEM_SOURCE_DIR}/reactions"
                      "${CHEM_SOURCE_DIR}/solvers"
                      "${FUNCS_SOURCE_DIR}")
set_property(GLOBAL PROPERTY CHEM_INCLUDES_DIR "${CHEM_INCLUDES_DIR}")
include_directories(${CHEM_INCLUDES_DIR})

set(GEOCHEM_LIBS geochemrxns geochembase geochemsolvers geochemutil)
set_property(GLOBAL PROPERTY GEOCHEM_LIBS "${GEOCHEM_LIBS}")

add_subdirectory(utility)
add_subdirectory(solvers)
add_subdirectory(base_chemistry)
add_subdirectory(reactions)

#
# Install targets
#
file(GLOB chem_inc "${CHEM_INCLUDES_DIR}/*")
add_install_include_file(${chem_inc})


if (BUILD_TESTS)
  get_property(CHEM_INCLUDES_DIR GLOBAL PROPERTY CHEM_INCLUDES_DIR)
  get_property(GEOCHEM_LIBS GLOBAL PROPERTY GEOCHEM_LIBS)

  include_directories(${CHEM_INCLUDES_DIR})
  include_directories(${UnitTest_INCLUDE_DIRS})

  file(GLOB chem_inc "${CHEM_INCLUDES_DIR}/*")
  add_install_include_file(${chem_inc})

  if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
    file(GLOB DataFiles "${CHEM_SOURCE_DIR}/test/*.xml"
                        "${CHEM_SOURCE_DIR}/test/*.dat"
                        "${CHEM_SOURCE_DIR}/test/native")
    file(COPY ${DataFiles} DESTINATION ${CHEM_BINARY_DIR}/test/)
  endif()

  add_amanzi_test(chemistry_native chemistry_native
                  KIND unit
                  SOURCE test/Main.cc   
                         test/chemistry_utilities.cc
                         test/chemistry_species.cc
                         test/chemistry_surface_site.cc           
                         test/chemistry_activity_model_factory.cc
                         test/chemistry_activity_coefficients.cc
                         test/chemistry_pitzer.cc
                         test/chemistry_secondary_species.cc
                         test/chemistry_aqueous_equilibrium_complex.cc
                         test/chemistry_mineral.cc
                         test/chemistry_mineral_kinetics_factory.cc
                         test/chemistry_kinetic_rate.cc
                         test/chemistry_beaker.cc           
                  LINK_LIBS
                         error_handling atk ${GEOCHEM_LIBS}
                         ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${UnitTest_LIBRARIES})

  add_amanzi_test(batch_native batch_native
                  KIND unit
                  SOURCE test/Main.cc
                         test/batch_native.cc
                  LINK_LIBS
                         error_handling atk ${GEOCHEM_LIBS}
                         ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${UnitTest_LIBRARIES})

  #
  # OpenMP test program
  #
  if (${OPENMP_FOUND})
    add_executable(beaker-threads test/beaker-threads.cc)
    target_link_libraries(beaker-threads ${GEOCHEM_LIBS} error_handling ${OpenMP_CXX_FLAGS})
    set_target_properties(beaker-threads PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
  endif()
endif()


