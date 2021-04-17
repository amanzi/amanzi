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
set(geochemrxns_src_files aqueous_equilibrium_complex.cc 
	                  beaker.cc beaker_constraints.cc
	                  general_rxn.cc 
	                  radioactive_decay.cc 
	                  ion_exchange_complex.cc 
	                  ion_exchange_rxn.cc 
	                  kinetic_rate.cc 
	                  kinetic_rate_tst.cc 
	                  Mineral.cc 
	                  mineral_kinetics_factory.cc 
	                  secondary_species.cc 
	                  simple_thermo_database.cc
	                  surface_complex.cc 
	                  surface_complexation_rxn.cc 
	                  sorption_isotherm.cc 
	                  sorption_isotherm_factory.cc 
	                  sorption_isotherm_freundlich.cc 
	                  sorption_isotherm_langmuir.cc 
	                  sorption_isotherm_linear.cc 
	                  sorption_isotherm_rxn.cc)

set(geochemrxns_inc_files aqueous_equilibrium_complex.hh
	                  beaker.hh
	                  general_rxn.hh
	                  radioactive_decay.hh
	                  ion_exchange_complex.hh
	                  ion_exchange_rxn.hh
	                  kinetic_rate.hh
	                  kinetic_rate_tst.hh
	                  Mineral.hh
	                  mineral_kinetics_factory.hh
	                  secondary_species.hh
	                  simple_thermo_database.hh
	                  surface_complex.hh
	                  surface_complexation_rxn.hh 
	                  sorption_isotherm.hh
	                  sorption_isotherm_factory.hh
	                  sorption_isotherm_freundlich.hh 
	                  sorption_isotherm_langmuir.hh
	                  sorption_isotherm_linear.hh
	                  sorption_isotherm_rxn.hh)

add_amanzi_library(geochemrxns SOURCE ${geochemrxns_src_files}
                               LINK_LIBS geochembase geochemsolvers geochemutil atk)

add_install_include_file(${geochemrxns_inc_files})
