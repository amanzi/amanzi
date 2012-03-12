#  -*- mode: cmake -*-

#
# Build TPL: Boost with CMake build 
#    
define_external_project(BoostCmake TARGET boost)

# ########################################################################### #
# Define Boost CMake build arguments                                          #
# ########################################################################### #

set(BoostCmake_CMAKE_ARGS)

# ########################################################################### #
# Add the Boost build target
# ########################################################################### #

# Make target: build boost with 'make boost' command
set(BoostCmake_target boost)

# Add the external project and tie to boost target
print_variable(BoostCmake_CMAKE_ARGS)
ExternalProject_Add(${BoostCmake_target}
    ${BoostCmake_ep_directory_args}
    ${BoostCmake_url_args}
    LIST_SEPARATOR ,
    CMAKE_ARGS
        ${Amanzi_CMAKE_COMPILER_ARGS}
	-DBUILD_PROJECTS:STRING=system,filesystem,program_options,regex
	-DBUILD_TESTS=OFF
	-DBUILD_TOOLS=OFF
	-DENABLE_SHARED=OFF
	-DENABLE_RELEASE=ON
	-DINSTALL_VERSIONED=OFF
	-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    ${BoostCmake_logging_opts}
)
