# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Base
#


set(SRCS 
	./activity_model_factory.cc
	./activity_model.cc
	./activity_model_debye_huckel.cc
	./activity_model_unit.cc
	./ion_exchange_site.cc
	./species.cc
	./surface_site.cc
)

add_library(geochembase ${SRCS}) 
add_install_library(geochembase)
link_list_add(geochembase)
