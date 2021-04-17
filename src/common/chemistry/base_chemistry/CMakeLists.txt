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
	                  VirialCoefficient.cc
	                  IonExchangeSite.cc
	                  Species.cc
	                  SurfaceSite.cc)

set(geochembase_inc_files ActivityModelFactory.hh
	                  ActivityModel.hh
	                  ActivityModelDebyeHuckel.hh
	                  ActivityModelPitzerHWM.hh
	                  ActivityModelUnit.hh
	                  VirialCoefficient.hh
	                  IonExchangeSite.hh
	                  Species.hh
	                  SurfaceSite.hh)

add_amanzi_library(geochembase SOURCE ${geochembase_src_files}
	                       LINK_LIBS geochemutil) 

add_install_include_file(${geochembase_inc_files})
