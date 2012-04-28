# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry_PK
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(TestManager)
include(InstallManager)

# Chemistry includes
include_directories(${CHEM_INCLUDES_DIR})

# Add the UnitTest includes
include_directories(${UnitTest_INCLUDE_DIRS})

#
# Install targets
#
file(GLOB chem_inc "${CHEM_INCLUDES_DIR}/*")
add_install_include_file(${chem_inc})

# Amanzi includes
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MESH_SIMPLE_SOURCE_DIR})
include_directories(${CHEMPK_SOURCE_DIR})

# Trilinos related (TPL) includes
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})

add_executable(chemistry_pk_tests
               main.cc
               test_chemistry_pk.cc)

target_link_libraries(chemistry_pk_tests chemistry_pk state 
                      simple_mesh error_handling exodus_reader  
                      ${UnitTest_LIBRARIES} 
                      ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES})
add_amanzi_test(chemistry_pk chemistry_pk_tests KIND unit)

if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
  execute_process(COMMAND ${CMAKE_COMMAND} -E
                  copy_directory ${CHEMPK_SOURCE_DIR}/tests_unit ${CHEMPK_BINARY_DIR}/tests_unit)
endif()


