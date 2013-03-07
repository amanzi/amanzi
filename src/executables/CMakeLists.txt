# -*- mode: cmake -*-

#
#  ATS
#    Executable
#

MACRO(WHOLE_ARCHIVE inlib)
  SET(${inlib} -Wl,--whole-archive ${inlib} -Wl,--no-whole-archive)
ENDMACRO(WHOLE_ARCHIVE inlib)

include_directories(${ATS_SOURCE_DIR}/src/data_structures)
include_directories(${ATS_SOURCE_DIR}/src/state)
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/mpc)

include_directories(${ATS_SOURCE_DIR}/src/coordinator)

add_executable(ats amanzi_unstructured_grid_simulation_driver.cc main.cc )

# wrap stuff in whole-archive for static stuff
WHOLE_ARCHIVE(pk_bases)
WHOLE_ARCHIVE(pk_test_divgrad)
WHOLE_ARCHIVE(pk_energy_advection_diffusion)
WHOLE_ARCHIVE(pk_energy_constant_temperature)
WHOLE_ARCHIVE(pk_energy_base)
WHOLE_ARCHIVE(pk_energy_three_phase)
WHOLE_ARCHIVE(pk_energy_two_phase)
WHOLE_ARCHIVE(pk_energy_surface_ice)
WHOLE_ARCHIVE(pk_transport_passive_tracer)
WHOLE_ARCHIVE(pk_flow_richards)
WHOLE_ARCHIVE(pk_flow_permafrost)
WHOLE_ARCHIVE(pk_flow_overland)
WHOLE_ARCHIVE(pk_flow_overland_head)
WHOLE_ARCHIVE(pk_flow_icy_overland)
WHOLE_ARCHIVE(pk_prescribed_deformation)
WHOLE_ARCHIVE(energy_relations_thermal_conductivity)
WHOLE_ARCHIVE(mpc)

WHOLE_ARCHIVE(relations_eos)
WHOLE_ARCHIVE(flow_relations_wrm)
WHOLE_ARCHIVE(flow_relations_wrm_models)
WHOLE_ARCHIVE(flow_relations_overland_conductivity)
WHOLE_ARCHIVE(flow_relations_surface_subsurface_fluxes)
WHOLE_ARCHIVE(flow_relations_elevation)
WHOLE_ARCHIVE(energy_relations_enthalpy)
WHOLE_ARCHIVE(energy_relations_internal_energy)
WHOLE_ARCHIVE(energy_relations_thermal_conductivity)
WHOLE_ARCHIVE(amanzi_mesh_functions)
WHOLE_ARCHIVE(amanzi_state)

target_link_libraries(ats
                      coordinator
                      ${mpc}
                      ${pk_test_divgrad}
#                      ${pk_transport_passive_tracer}
                      ${pk_energy_constant_temperature}
                      ${pk_energy_advection_diffusion}
                      ${pk_energy_base}
                      ${pk_energy_two_phase}
                      ${pk_energy_three_phase}
                      ${pk_energy_surface_ice}
                      ${pk_flow_permafrost}
                      ${pk_flow_richards}
                      ${pk_flow_overland}
                      ${pk_flow_overland_head}
                      ${pk_flow_icy_overland}
		              ${pk_prescribed_deformation}
                      ${pk_bases}
                      ${energy_relations_thermal_conductivity}
                      ${energy_relations_internal_energy}
                      ${energy_relations_enthalpy}
                      ${flow_relations_wrm}
                      ${flow_relations_wrm_models}
                      ${flow_relations_overland_conductivity}
                      ${flow_relations_elevation}
                      ${flow_relations_surface_subsurface_fluxes}
                      ${relations_eos}
                      bc_factory
                      advection
                      divgrad
                      amanzi_whetstone
                      amanzi_solvers
                      amanzi_time_integration
                      amanzi_mesh
                      amanzi_mesh_audit
                      amanzi_mesh_functions
                      amanzi_geometry
                      amanzi_mesh_factory
                      amanzi_input
                      amanzi_interface
                      amanzi_output
                      amanzi_stk_mesh
                      ${amanzi_state}
                      )

install(TARGETS ats DESTINATION bin)






