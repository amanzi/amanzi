include_directories(${PROJECT_SOURCE_DIR}/includes)

#
# library
#
set(SRCS 
	./AqueousEquilibriumComplex.cpp 
	./Beaker.cpp 
	./GeneralRxn.cpp 
	./IonExchangeComplex.cpp 
	./KineticRate.cpp 
	./KineticRateTST.cpp 
	./Mineral.cpp 
	./MineralKineticsFactory.cpp 
	./Reaction.cpp 
	./SecondarySpecies.cpp 
	./SimpleThermoDatabase.cpp
	./SurfaceComplex.cpp 
	./SurfaceComplexationRxn.cpp 
)

add_library(geochemrxns ${SRCS}) 

#
# glenn's test program
#
set(EXE_NAME simpleMain)

set(EXE_SRCS
	./SimpleMain.cpp
)

set(GEOCHEM_LIBS
	geochemrxns
	geochembase
	geochemsolvers
	geochemutil
)

add_executable(${EXE_NAME} ${EXE_SRCS})
target_link_libraries(${EXE_NAME} ${GEOCHEM_LIBS})
