# -*- mode: cmake -*-

#
#  Amanzi
#    Chemistry Base
#


set(SRCS 
	./ActivityModelFactory.cpp 
	./ActivityModel.cpp 
	./ActivityModelDebyeHuckel.cpp 
	./ActivityModelUnit.cpp 
	./IonExchangeSite.cpp 
	./Species.cpp 
	./SurfaceSite.cpp 
)

add_library(geochembase ${SRCS}) 

