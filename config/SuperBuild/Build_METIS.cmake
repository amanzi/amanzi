#  -*- mode: cmake -*-

#
# Build METIS 
#    

include(ExternalProject)
include(BuildLibraryName)

include(TPLVersions)


# Define source, build and install directories
set(METIS_source_dir "${CMAKE_BINARY_DIR}/external-projects/metis/metis-${METIS_VERSION}-source")
set(METIS_binary_dir "${CMAKE_BINARY_DIR}/external-projects/metis/metis-${METIS_VERSION}-build")
set(METIS_install_dir "${CMAKE_BINARY_DIR}/external-projects/metis")

# Make target: build metis with 'make metis' command
set(METIS_target metis)

# Add the external project and tie to metis target
set(METIS_GKlib_path ${METIS_source_dir}/GKlib)
ExternalProject_Add(${METIS_target}
    SOURCE_DIR ${METIS_source_dir}
    BINARY_DIR ${METIS_binary_dir}
    INSTALL_DIR ${METIS_install_dir}
    URL ${METIS_URL_STRING}/${METIS_ARCHIVE_FILE}
    URL_MD5 ${METIS_MD5_SUM}
    CMAKE_ARGS
            ${Amanzi_CMAKE_COMPILER_ARGS}
	    -D CMAKE_INSTALL_PREFIX:FILEPATH=${METIS_install_dir}
	    -D GKLIB_PATH:FILEPATH=${METIS_GKlib_path}
	    -D SHARED:Bool=${BUILD_SHARED_LIBS}
)
                                      
# Define variables needed by other packages
# These should match the final output from FindMETIS.cmake

# Include directories
set(METIS_INCLUDE_DIR "${METIS_install_dir}/include")
set(METIS_INCLUDE_DIRS "${METIS_install_dir}/include")

# Libraries
build_library_name(metis METIS_LIBRARY_FILENAME)
set(METIS_LIBRARY ${METIS_install_dir}/lib/${METIS_LIBRARY_FILENAME})
set(METIS_LIBRARIES "${METIS_LIBRARY}")
