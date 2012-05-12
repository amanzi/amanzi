# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Test Drivers
#

# Chemistry includes
include_directories(${CHEM_INCLUDES_DIR})

#
# batch chemistry drivers and test program
#

set(BATCH_EXE_NAME batch_chem)

set(BATCH_EXE_SRCS batch_chem.cc)

set(GEOCHEM_LIBS
	geochemrxns
	geochembase
	geochemsolvers
	geochemutil
)


add_executable(${BATCH_EXE_NAME} ${BATCH_EXE_SRCS})
target_link_libraries(${BATCH_EXE_NAME} ${GEOCHEM_LIBS} error_handling)

if (NOT (${CHEM_SOURCE_DIR} EQUAL ${CHEM_BINARY_DIR}))
  execute_process(COMMAND ${CMAKE_COMMAND} -E 
                  copy_directory ${CHEM_SOURCE_DIR}/tests_regression/input 
                  ${CHEM_BINARY_DIR}/tests_regression/input) 
endif()


#
# OpenMP test program
#
if (${OPENMP_FOUND})
  add_executable(beaker-threads beaker-threads.cc)
  target_link_libraries(beaker-threads ${GEOCHEM_LIBS} error_handling ${OpenMP_CXX_FLAGS})
  set_target_properties(beaker-threads PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
endif()

