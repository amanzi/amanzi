# -*- mode: cmake -*-

#
#  ATS
#    Executable
#

include_directories(${ATS_SOURCE_DIR}/src/data_structures)
include_directories(${ATS_SOURCE_DIR}/src/state)
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/mpc)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad/upwind_scheme)
include_directories(${ATS_SOURCE_DIR}/src/operators/advection)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/base)


include_directories(${ATS_SOURCE_DIR}/src/coordinator)

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
include_evaluators_directories(LISTNAME TEST_PKS_DIVGRAD_TEST_REG_INCLUDES)
include_evaluators_directories(LISTNAME MPC_PK_FAC_INCLUDES)
#include_evaluators_directories(LISTNAME TRANSPORT_PASSIVE_TRACER_REG_INCLUDES)


add_executable(ats amanzi_unstructured_grid_simulation_driver.cc main.cc )


target_link_libraries(ats
                      coordinator
                      mpc
                      pk_test_divgrad
#                      pk_transport_passive_tracer
                      pk_energy_constant_temperature
                      pk_energy_advection_diffusion
		              pk_energy_two_phase
                      pk_energy_three_phase
                      pk_energy_surface_ice
                      pk_energy_base
		              pk_flow_permafrost
                      pk_flow_richards
                      pk_flow_overland
                      pk_flow_overland_head
                      pk_flow_icy_overland
		              pk_prescribed_deformation
#                      pk_prescribed_volumetric_deformation
                      pk_volumetric_deformation
                      pk_surface_balance_SEB
                      pk_bases
                      deform_relations_porosity
                      energy_relations_thermal_conductivity
                      energy_relations_internal_energy
                      energy_relations_enthalpy
                      energy_relations_source_terms
                      flow_relations_wrm
                      flow_relations_wrm_models
                      flow_relations_overland_conductivity
                      flow_relations_elevation
                      flow_relations_porosity
                      flow_relations_surface_subsurface_fluxes
                      relations_eos
                      relations_ewc
                      bc_factory
                      advection
                      divgrad
                      deformation_operator
                      amanzi_state
		              amanzi_whetstone
                      amanzi_time_integration
		              amanzi_solvers
                      amanzi_mesh
                      amanzi_mesh_audit
                      amanzi_mesh_functions
                      amanzi_geometry
                      amanzi_mesh_factory
                      amanzi_input
                      amanzi_interface
                      amanzi_output
                      amanzi_stk_mesh
                      )

install(TARGETS ats DESTINATION bin)
