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
# Included souce
#------------------------------------------------------------------------------#

include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/mpc)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad/upwind_scheme)
include_directories(${ATS_SOURCE_DIR}/src/operators/advection)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/base)
include_directories(${ATS_SOURCE_DIR}/src/coordinator)

#------------------------------------------------------------------------------#
# ATS generated include files
#------------------------------------------------------------------------------#

include_evaluators_directories(LISTNAME ATS_POROSITY_EVAL_INCLUDES)
include_evaluators_directories(LISTNAME CONSTITUTIVE_RELATIONS_EOS_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME CONSTITUTIVE_RELATIONS_SURFACE_SUBSURFACE_FLUXES_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME FLOW_CONSTITUTIVE_RELATIONS_OVERLAND_CONDUCTIVITY_INCLUDES)
include_evaluators_directories(LISTNAME FLOW_CONSTITUTIVE_RELATIONS_POROSITY_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME WRM_MODELS_EVALUATORS_FACTORIES_INCLUDES)
include_evaluators_directories(LISTNAME WMR_EVALUATORS_FACTORIES_INCLUDES)
include_evaluators_directories(LISTNAME FLOW_ICY_OVERLAND_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME FLOW_OVERLAND_HEAD_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME FLOW_OVERLAND_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME RICHARDS_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME PERMAFROST_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME PRESCRIBED_DEFORMATION_EVAL_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_ADVECTION_DIFFUSION_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_CONSTANT_TEMPERATURE_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_CONST_REL_INTERNAL_ENERGY_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_CONST_REL_SOURCE_TERMS_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_CONST_REL_THERMAL_CONDUCTIVITY_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_SURFACE_ICE_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_TWO_PHASE_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_THREE_PHASE_REG_INCLUDES)
include_evaluators_directories(LISTNAME SURFACE_BALANCE_SEB_REG_INCLUDES)
include_evaluators_directories(LISTNAME TEST_PKS_REG_INCLUDES)
include_evaluators_directories(LISTNAME MPC_PK_FAC_INCLUDES)
#include_evaluators_directories(LISTNAME TRANSPORT_PASSIVE_TRACER_REG_INCLUDES)


#------------------------------------------------------------------------------#
# ATS cross language library
#------------------------------------------------------------------------------#

add_library(
	# Library name
	atsx

	# C source
	c/ats_interface.cc
	c/ats_clm_driver.cc

	# Fortran source
	c/ats_f90.c
	fc/ats_data.f90
	fc/ats_defines.f90
	fc/ats_bindings.f90
	fc/ats_interface.f90
	fc/ats.f90
)

install(TARGETS atsx DESTINATION lib)
#install(FILES ats.mod DESTINATION lib)

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

