include_directories(${PROJECT_SOURCE_DIR}/includes)
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
target_link_libraries(${EXE_NAME} ${GEOCHEM_LIBS})

if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
  execute_process(COMMAND ${CMAKE_COMMAND} -E 
                  copy_directory ${CMAKE_SOURCE_DIR}/chemistry/test_drivers/input 
                  ${CMAKE_BINARY_DIR}/chemistry/test_drivers/input) 
endif()


#INSTALL(TARGETS ${EXE_NAME} DESTINATION ${CMAKE_HOME_DIRECTORY}/bin)
