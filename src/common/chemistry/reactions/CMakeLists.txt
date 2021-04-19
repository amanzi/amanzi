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
set(geochemrxns_src_files AqueousEquilibriumComplex.cc 
	                  Beaker.cc Beaker_Constraints.cc
	                  GeneralRxn.cc 
	                  RadioactiveDecay.cc 
	                  IonExchangeComplex.cc 
	                  IonExchangeRxn.cc 
	                  KineticRate.cc 
	                  KineticRateFactory.cc 
	                  KineticRateTST.cc 
	                  Mineral.cc 
	                  SecondarySpecies.cc 
	                  SimpleThermoDatabase.cc
	                  SurfaceComplex.cc 
	                  SurfaceComplexationRxn.cc 
	                  SorptionIsotherm.cc 
	                  SorptionIsothermFactory.cc 
	                  SorptionIsothermFreundlich.cc 
	                  SorptionIsothermLangmuir.cc 
	                  SorptionIsothermLinear.cc 
	                  SorptionIsothermRxn.cc)

set(geochemrxns_inc_files AqueousEquilibriumComplex.hh
	                  Beaker.hh
	                  BeakerState.hh
	                  GeneralRxn.hh
	                  RadioactiveDecay.hh
	                  IonExchangeComplex.hh
	                  IonExchangeRxn.hh
	                  KineticRate.hh
	                  KineticRateFactory.hh
	                  KineticRateTST.hh
	                  Mineral.hh
	                  SecondarySpecies.hh
	                  SimpleThermoDatabase.hh
	                  SurfaceComplex.hh
	                  SurfaceComplexationRxn.hh 
	                  SorptionIsotherm.hh
	                  SorptionIsothermFactory.hh
	                  SorptionIsothermFreundlich.hh 
	                  SorptionIsothermLangmuir.hh
	                  SorptionIsothermLinear.hh
	                  SorptionIsothermRxn.hh)

add_amanzi_library(geochemrxns SOURCE ${geochemrxns_src_files}
                               LINK_LIBS geochembase geochemsolvers geochemutil
                                         atk functions)

add_install_include_file(${geochemrxns_inc_files})
