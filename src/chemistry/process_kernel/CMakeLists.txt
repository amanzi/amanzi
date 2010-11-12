include_directories(${MPI_INCLUDES})
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/mpc)
include_directories(${PROJECT_SOURCE_DIR}/simple_mesh)
include_directories(${PROJECT_SOURCE_DIR}/chemistry/includes)

set(GEOCHEM_LIBS
	geochemrxns
	geochembase
	geochemsolvers
	geochemutil
)

add_library(chemistry Chemistry_PK.cpp)

target_link_libraries(chemistry ${GEOCHEM_LIBS} ${Trilinos_LIBRARIES})

link_directories(${Trilinos_LIBRARY_DIRS})




