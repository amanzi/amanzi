# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Test Drivers
#

# Chemistry includes
include_directories(${CHEM_INCLUDES_DIR})

#
# ben's test program
#
set(EXE_NAME geochemistryTest)

set(EXE_SRCS 
	./geochemistry_test.cc
)

set(GEOCHEM_LIBS
	geochemrxns
	geochembase
	geochemsolvers
	geochemutil
)


add_executable(${EXE_NAME} ${EXE_SRCS})
target_link_libraries(${EXE_NAME} ${GEOCHEM_LIBS} error_handling)

if (NOT (${CHEM_SOURCE_DIR} EQUAL ${CHEM_BINARY_DIR}))
  execute_process(COMMAND ${CMAKE_COMMAND} -E 
                  copy_directory ${CHEMPK_SOURCE_DIR}/tests_regression/input 
                  ${CHEMPK_BINARY_DIR}/tests_regression/input) 
endif()


#
# OpenMP test program
#
if (${OPENMP_FOUND})
  add_executable(beaker-threads beaker-threads.cc)
  target_link_libraries(beaker-threads ${GEOCHEM_LIBS} error_handling ${OpenMP_CXX_FLAGS})
  set_target_properties(beaker-threads PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
endif()

