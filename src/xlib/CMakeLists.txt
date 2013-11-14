#------------------------------------------------------------------------------#
# ATS
#
# License: see $ATS_DIR/COPYRIGHT
#
# CMake rules for C and Fortran Language interfaces.
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------#

include_directories(${ATS_SOURCE_DIR}/src/xlib/c)
include_directories(${ATS_SOURCE_DIR}/src/xlib/fc)

#------------------------------------------------------------------------------#
# ATS cross language library
#------------------------------------------------------------------------------#

add_library(
	# Library name
	atsx

	# C source
	c/ats_interface.cc

	# Fortran source
	c/ats_f90.c
	fc/ats_data.f90
	fc/ats_defines.f90
	fc/ats_bindings.f90
	fc/ats_interface.f90
	fc/ats.f90
)

#------------------------------------------------------------------------------#
# Command to update ATS defines
#------------------------------------------------------------------------------#

add_custom_command(
	TARGET atsx
	PRE_BUILD
	COMMAND ${ATS_SOURCE_DIR}/tools/make_defines ${ATS_SOURCE_DIR}/src/xlib
	DEPENDS ${ATS_SOURCE_DIR}/src/xlib/ats_defines
)

#------------------------------------------------------------------------------#
# Don't link MPI CXX libraries.
#------------------------------------------------------------------------------#

set_property(TARGET atsx PROPERTY
	COMPILE_DEFINITIONS OMPI_SKIP_MPICXX MPICH_SKIP_MPICXX
)
