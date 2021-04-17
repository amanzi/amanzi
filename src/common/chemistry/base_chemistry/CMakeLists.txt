# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Base
#

set(geochembase_src_files ActivityModelFactory.cc
	                  ActivityModel.cc
	                  ActivityModelDebyeHuckel.cc
	                  ActivityModelPitzerHWM.cc 
	                  ActivityModelUnit.cc
	                  virial_coefficient.cc
	                  ion_exchange_site.cc
	                  species.cc
	                  surface_site.cc)

set(geochembase_inc_files ActivityModelFactory.hh
	                  ActivityModel.hh
	                  ActivityModelDebyeHuckel.hh
	                  ActivityModelPitzerHWM.hh
	                  ActivityModelUnit.hh
	                  virial_coefficient.hh
	                  ion_exchange_site.hh
	                  species.hh
	                  surface_site.hh)

add_amanzi_library(geochembase SOURCE ${geochembase_src_files}
	                       LINK_LIBS geochemutil) 

add_install_include_file(${geochembase_inc_files})
