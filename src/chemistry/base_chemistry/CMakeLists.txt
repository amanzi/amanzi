include_directories(${PROJECT_SOURCE_DIR}/includes)

set(SRCS 
	./ActivityModelFactory.cpp 
	./ActivityModel.cpp 
	./ActivityModelDebyeHuckel.cpp 
	./ActivityModelUnit.cpp 
	./IonExchangeSite.cpp 
	./Species.cpp 
)

add_library(geochembase ${SRCS}) 

