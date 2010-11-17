include_directories(${PROJECT_SOURCE_DIR}/includes)

set(SRCS 
	 ./Block.cpp 
	 ./FileIO.cpp 
	 ./Glenn.cpp 
	 ./ChemistryException.cpp
	 ./StringTokenizer.cpp
)

add_library(geochemutil ${SRCS}) 

