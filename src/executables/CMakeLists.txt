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
WHOLE_ARCHIVE(pk_factory)
WHOLE_ARCHIVE(pk_energy_two_phase)
WHOLE_ARCHIVE(pk_energy_advection_diffusion)
WHOLE_ARCHIVE(pk_energy_constant_temperature)
WHOLE_ARCHIVE(pk_energy_two_phase)
WHOLE_ARCHIVE(pk_transport_passive_tracer)
WHOLE_ARCHIVE(pk_flow_richards)
WHOLE_ARCHIVE(energy_relations_thermal_conductivity)
WHOLE_ARCHIVE(mpc)

WHOLE_ARCHIVE(flow_relations_eos)
WHOLE_ARCHIVE(flow_relations_wrm)

target_link_libraries(ats
                      coordinator
                      ${mpc}
                      ${pk_factory}
                      ${pk_transport_passive_tracer}
                      ${pk_energy_constant_temperature}
                      ${pk_energy_advection_diffusion}
                      ${pk_energy_two_phase}
                      ${pk_flow_richards}
                      energy_bc_factory
                      ${energy_relations_thermal_conductivity}
                      energy_relations_internal_energy
                      flow_bc_factory
                      ${flow_relations_eos}
                      ${flow_relations_wrm}
                      advection
                      divgrad
                      amanzi_whetstone
                      amanzi_solvers
                      amanzi_time_integration
                      amanzi_mesh
                      amanzi_mesh_functions
                      amanzi_cell_geometry
                      amanzi_mesh_factory
                      amanzi_input
                      amanzi_interface
                      amanzi_output
                      amanzi_stk_mesh
                      amanzi_state
                      )

install(TARGETS ats DESTINATION ${ATS_SOURCE_DIR}/bin)






