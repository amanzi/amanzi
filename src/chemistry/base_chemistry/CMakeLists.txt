# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Base
#


set(SRCS 
	./activity-model-factory.cc
	./activity-model.cc
	./activity-model-debye-huckel.cc
	./activity-model-unit.cc
	./IonExchangeSite.cpp
	./Species.cpp 
	./SurfaceSite.cpp 
)

add_library(geochembase ${SRCS}) 
add_install_library(geochembase)

