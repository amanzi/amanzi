# -*- mode: cmake -*-
# Amanzi Version Information
#
include(MercurialMacros)
include(PrintVariable)
include(InstallManager)


# Amanzi version format major.minor.patch
set(AMANZI_VERSION_MAJOR 0)
set(AMANZI_VERSION_MINOR 1)
set(AMANZI_VERSION_PATCH 0)
set(AMANZI_VERSION "0.1")

# Set Mercurial branch and ids
mercurial_branch(AMANZI_HG_BRANCH)
mercurial_global_id(AMANZI_HG_GLOBAL_HASH)
mercurial_local_id(AMANZI_HG_LOCAL_ID)


# Status output
message(STATUS "Amanzi Version:\t${AMANZI_VERSION}")
message(STATUS "\tMercurial Branch:\t${AMANZI_HG_BRANCH}")
message(STATUS "\tMercurial Global ID:\t${AMANZI_HG_GLOBAL_HASH}")
message(STATUS "\tMercurial Local ID:\t${AMANZI_HG_LOCAL_ID}")

# Write the version header file
set(version_template ${AMANZI_SOURCE_TOOLS_DIR}/cmake/amanzi_version.hh.in)
configure_file(${version_template}
               ${CMAKE_CURRENT_BINARY_DIR}/amanzi_version.hh
               @ONLY)

add_install_include_file(${CMAKE_CURRENT_BINARY_DIR}/amanzi_version.hh)             

             


