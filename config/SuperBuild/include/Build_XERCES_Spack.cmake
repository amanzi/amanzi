#  -*- mode: cmake -*-
#
# Build TPL: XERCES using Spack
#  

# -- Define all the directories and common external project flags
define_external_project_args(XERCES TARGET xerces BUILD_IN_SOURCE)

ExternalProject_Add(${XERCES_BUILD_TARGET} # Xerces
    
    # The following 2 lines are only here because CMake will not allow
    # a build without at least the first three variables defined.
    # However, they serve no real purpose for this particular build.
    DOWNLOAD_COMMAND ls
    
    CONFIGURE_COMMAND ls
    
    BUILD_COMMAND ${SPACK_BINARY} install xerces-c 
    
    INSTALL_COMMAND  ${SPACK_BINARY} view symlink ${TPL_INSTALL_PREFIX} xerces-c 
)

build_library_name(xerces-c XERCES_LIBRARY APPEND_PATH ${TPL_INSTALL_PREFIX}/lib)

set(XERCES_LIBRARIES ${XERCES_LIBRARY})
set(XERCES_INCLUDE_DIRS ${TPL_INSTALL_PREFIX}/include)
set(XERCES_INSTALL_PREFIX ${TPL_INSTALL_PREFIX})
append_set(Amanzi_TPL_CMAKE_ARGS
    -DXERCES_DIR:FILEPATH=${XERCES_INSTALL_PREFIX})

amanzi_tpl_version_write(FILENAME ${TPL_VERSIONS_INCLUDE_FILE}
    PREFIX XERCES
    VERSION ${TMP_SPACK_VERSION} )

