#  -*- mode: cmake -*-

#
# Build TPL:  CGNS 
#    
define_external_project(CGNS TARGET cgns 
                             DEPENDS curl
                             BUILD_IN_SOURCE)

# Add the external project and tie to zlib target
ExternalProject_Add(${CGNS_target}
    DEPENDS curl
    ${CGNS_ep_directory_args}
    ${CGNS_url_args}
    CONFIGURE_COMMAND
                    <SOURCE_DIR>/configure 
                              --prefix=<INSTALL_DIR>
                              --enable-lfs
                              --enable-64bit
                              --with-fortran=no
                              FC=false
    BUILD_IN_SOURCE TRUE      
    ${CGNS_logging_opts}               
)
