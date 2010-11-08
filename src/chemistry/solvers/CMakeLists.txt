include_directories(${PROJECT_SOURCE_DIR}/includes)

set(SRCS 
	 ./LU.cpp 
	 ./Newton.cpp 
)

add_library(geochemsolvers ${SRCS}) 

