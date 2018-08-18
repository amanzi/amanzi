# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Base
#

set(geochembase_src_files activity_model_factory.cc
	                  activity_model.cc
	                  activity_model_debye_huckel.cc
	                  activity_model_pitzer_hwm.cc 
	                  activity_model_unit.cc
	                  virial_coefficient.cc
	                  ion_exchange_site.cc
	                  species.cc
	                  surface_site.cc)

set(geochembase_inc_files activity_model_factory.hh
	                  activity_model.hh
	                  activity_model_debye_huckel.hh
	                  activity_model_pitzer_hwm.hh
	                  activity_model_unit.hh
	                  virial_coefficient.hh
	                  ion_exchange_site.hh
	                  species.hh
	                  surface_site.hh)

add_amanzi_library(geochembase SOURCE ${geochembase_src_files}
	                       LINK_LIBS geochemutil) 

add_install_include_file(${geochembase_inc_files})
