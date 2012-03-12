#  -*- mode: cmake -*-

#
# Build TPL:  CCSE 
#    
define_external_project(CCSE TARGET ccse)

# With ENABLE_MPI CCSE requires MPI wrappers and sets compilers to
# <MPI_PREFIX>/bin/mpi* in the root CMakeLists.txt
#if ( NOT MPI_PREFIX )
#  message(FATAL_ERROR "When compiling CCSE you must define an root MPI installation"
#                      " path with the variable MPI_PREFIX."
#		      "\n-DMPI_PREFIX:FILEPATH=<MPI Root installation path>\n"#
#"Please define this argument and re-run cmake."
#  )
#endif()

# ############################################################################ #
# Comman CMake arguments for 2D and 3D
# ############################################################################ #
message(STATUS "Build CCSE with space dimension ${CCSE_BL_SPACEDIM}")
set(CCSE_CMAKE_COMMON_ARGS
                       -DENABLE_Config_Report:BOOL=TRUE
		       -DENABLE_MPI:BOOL=TRUE
		       -DMPI_PREFIX:FILEPATH=${MPI_PREFIX}
		       -DENABLE_OpenMP:BOOL=TRUE
		       -DENABLE_TESTS:BOOL=FALSE
		       -DBL_PRECISION:STRING=DOUBLE
		       -DBL_SPACEDIM:INT=${CCSE_BL_SPACEDIM})
		     
# ############################################################################ #
# Add external project
# ############################################################################ #
ExternalProject_Add(${CCSE_target}
    ${CCSE_ep_directory_args}
    ${CCSE_url_args}
    CMAKE_ARGS
             ${Amanzi_CMAKE_COMPILER_ARGS}
	     ${CCSE_CMAKE_COMMON_ARGS}
	     -DCMAKE_INSTALL_PREFIX:FILEPATH=${TPL_INSTALL_PREFIX}
    ${CCSE_logging_opts}                  
)
