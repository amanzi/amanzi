include_directories(${PROJECT_SOURCE_DIR}/includes)
#
# ben's test program
#
set(EXE_NAME geochemistryTest)

set(EXE_SRCS 
	./ThermoDatabase.cpp 
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

INSTALL(TARGETS ${EXE_NAME} DESTINATION ${CMAKE_HOME_DIRECTORY}/bin)
