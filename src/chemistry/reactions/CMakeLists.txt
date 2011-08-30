# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Reactions
#

# Chemistry includes
include_directories(${CHEM_INCLUDES_DIR})


#
# library
#
set(SRCS 
	./aqueous_equilibrium_complex.cc 
	./beaker.cc 
	./general_rxn.cc 
	./ion_exchange_complex.cc 
	./kinetic_rate.cc 
	./kinetic_rate_tst.cc 
	./mineral.cc 
	./mineral_kinetics_factory.cc 
	./secondary_species.cc 
	./simple_thermo_database.cc
	./surface_complex.cc 
	./surface_complexation_rxn.cc 
	./sorption_isotherm.cc 
	./sorption_isotherm_factory.cc 
	./sorption_isotherm_freundlich.cc 
	./sorption_isotherm_langmuir.cc 
	./sorption_isotherm_linear.cc 
	./sorption_isotherm_rxn.cc 
)

add_library(geochemrxns ${SRCS}) 
add_install_library(geochemrxns)

#
# glenn's test program
#
set(EXE_NAME simpleMain)

set(EXE_SRCS
	./simple_main.cc
)

set(GEOCHEM_LIBS
	geochemrxns
	geochembase
	geochemsolvers
	geochemutil
)

add_executable(${EXE_NAME} ${EXE_SRCS})
target_link_libraries(${EXE_NAME} ${GEOCHEM_LIBS} error_handling)
