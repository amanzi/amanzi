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
	./geochemistry_test.cpp 
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
                  copy_directory ${CHEM_SOURCE_DIR}/test_drivers/input 
                  ${CHEM_BINARY_DIR}/test_drivers/input) 
endif()


#
# openmp example
#
#add_executable(beaker-threads beaker-threads.cc)
#target_link_libraries(beaker-threads ${GEOCHEM_LIBS} error_handling gomp)
#set_target_properties(beaker-threads PROPERTIES COMPILE_FLAGS "-fopenmp")

